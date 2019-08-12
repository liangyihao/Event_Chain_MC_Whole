/***************************************************************************
 *cr                                                                       
 *cr            (C) Copyright 1995-2016 The Board of Trustees of the           
 *cr                        University of Illinois                       
 *cr                         All Rights Reserved                        
 *cr                                                                   
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: dcdplugin.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.79 $       $Date: 2016/11/28 05:01:53 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Code for reading and writing CHARMM, NAMD, and X-PLOR format 
 *   molecular dynamic trajectory files.
 *
 *  Standalone test binary compilation flags:
 *  cc -fast -xarch=v9a -I../../include -DTEST_DCDPLUGIN dcdplugin.c \
 *    -o ~/bin/readdcd -lm
 *
 *  Profiling flags:
 *  cc -xpg -fast -xarch=v9a -g -I../../include -DTEST_DCDPLUGIN dcdplugin.c \
 *    -o ~/bin/readdcd -lm
 *
 ***************************************************************************/

#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "dcd_writer.hpp"
#include "public.hpp"
#include "basic.hpp"
#include "manage.hpp"
#include "Input_File_Parser.hpp"
/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2016 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: largefiles.h,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.4 $       $Date: 2016/11/28 05:01:54 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Platform dependent defines for enabling 64-bit file I/O on 32-bit machines
 *
 ***************************************************************************/
 
#if defined(_AIX)
/* Define to enable large file extensions on AIX */
#define _LARGE_FILE
#define _LARGE_FILES
#else
/* Defines which enable LFS I/O interfaces for large (>2GB) file support
 * on 32-bit machines.  These must be defined before inclusion of any
 * system headers.
 */
#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE
#endif
#define _FILE_OFFSET_BITS 64
#endif



#include "dcd_fastio.hpp"       /* must come before others, for O_DIRECT...   */

/*
 * Per-timestep atom coordinates and periodic cell information
 */ 
typedef struct {
  float *coords;        /**< coordinates of all atoms, arranged xyzxyzxyz   */
  float *velocities;    /**< space for velocities of all atoms; same layout */
                        /**< NULL unless has_velocities is set              */

  /*@{*/   
  /**
   * Unit cell specification of the form A, B, C, alpha, beta, gamma.
   * notes: A, B, C are side lengths of the unit cell
   * alpha = angle between b and c
   *  beta = angle between a and c
   * gamma = angle between a and b
   */ 
  float A, B, C, alpha, beta, gamma; 
  /*@}*/
  double physical_time; /**< physical time point associated with this frame */
} molfile_timestep_t;
#define MOLFILE_ERROR            -1   /**< error reading/opening a file   */
#define MOLFILE_SUCCESS           0   /**< succeeded in reading file      */

/////////
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661922
#endif

#define RECSCALE32BIT 1
#define RECSCALE64BIT 2
#define RECSCALEMAX   2

typedef struct {
  fio_fd fd;
  int natoms;//Number of atoms per frame
  int nsets;//Number of frames
  int setsread;
  int istart;//Starting timestep of DCD file
  int nsavc;//Number of timesteps between written DCD frames
  double delta;//length of a timestep
  int nfixed;
  float *x, *y, *z;
  int *freeind;
  float *fixedcoords;
  int reverse;
  int charmm;  
  int first;
  int with_unitcell;
} dcdhandle;

/* Define error codes that may be returned by the DCD routines */
#define DCD_SUCCESS      0  /* No problems                     */
#define DCD_EOF         -1  /* Normal EOF                      */
#define DCD_DNE         -2  /* DCD file does not exist         */
#define DCD_OPENFAILED  -3  /* Open of DCD file failed         */
#define DCD_BADREAD     -4  /* read call on DCD file failed    */
#define DCD_BADEOF      -5  /* premature EOF found in DCD file */
#define DCD_BADFORMAT   -6  /* format of DCD file is wrong     */
#define DCD_FILEEXISTS  -7  /* output file already exists      */
#define DCD_BADMALLOC   -8  /* malloc failed                   */
#define DCD_BADWRITE    -9  /* write call on DCD file failed   */

/* Define feature flags for this DCD file */
#define DCD_IS_XPLOR        0x00
#define DCD_IS_CHARMM       0x01
#define DCD_HAS_4DIMS       0x02
#define DCD_HAS_EXTRA_BLOCK 0x04
#define DCD_HAS_64BIT_REC   0x08

/* defines used by write_dcdstep */
#define NFILE_POS 8L
#define NSTEP_POS 20L

/* READ Macro to make porting easier */
#define READ(fd, buf, size)  fio_fread(((void *) buf), (size), 1, (fd))

/* WRITE Macro to make porting easier */
#define WRITE(fd, buf, size) fio_fwrite(((void *) buf), (size), 1, (fd))

/* XXX This is broken - fread never returns -1 */
#define CHECK_FREAD(X, msg) if (X==-1) { return(DCD_BADREAD); }
#define CHECK_FEOF(X, msg)  if (X==0)  { return(DCD_BADEOF); }



/* print DCD error in a human readable way */
static void print_dcderror(const char *func, int errcode) {
  const char *errstr;

  switch (errcode) {
    case DCD_EOF:         errstr = "end of file"; break;
    case DCD_DNE:         errstr = "file not found"; break;
    case DCD_OPENFAILED:  errstr = "file open failed"; break;
    case DCD_BADREAD:     errstr = "error during read"; break;
    case DCD_BADEOF:      errstr = "premature end of file"; break;
    case DCD_BADFORMAT:   errstr = "corruption or unrecognized file structure"; break;
    case DCD_FILEEXISTS:  errstr = "output file already exists"; break;
    case DCD_BADMALLOC:   errstr = "memory allocation failed"; break;
    case DCD_BADWRITE:    errstr = "error during write"; break;
    case DCD_SUCCESS:     
    default:
      errstr = "no error";
      break;
  } 
  printf("dcdplugin) %s: %s\n", func, errstr); 
}

/* 
 * Write a timestep to a dcd file
 * Input: fd - a file struct for which a dcd header has already been written
 *       curframe: Count of frames written to this file, starting with 1.
 *       curstep: Count of timesteps elapsed = istart + curframe * nsavc.
 *        natoms - number of elements in x, y, z arrays
 *        x, y, z: pointers to atom coordinates
 * Output: 0 on success, negative error code on failure.
 * Side effects: coordinates are written to the dcd file.
 */
static int write_dcdstep(fio_fd fd, int curframe, int curstep, int N, 
                  const float *X, const float *Y, const float *Z, 
                  const double *unitcell, int charmm) {
  int out_integer;

  if (charmm) {
    /* write out optional unit cell */
    if (unitcell != NULL) {
      out_integer = 48; /* 48 bytes (6 floats) */
      fio_write_int32(fd, out_integer);
      WRITE(fd, unitcell, out_integer);
      fio_write_int32(fd, out_integer);
    }
  }

  /* write out coordinates */
  out_integer = N*4; /* N*4 bytes per X/Y/Z array (N floats per array) */
  fio_write_int32(fd, out_integer);
  if (fio_fwrite((void *) X, out_integer, 1, fd) != 1) return DCD_BADWRITE;
  fio_write_int32(fd, out_integer);
  fio_write_int32(fd, out_integer);
  if (fio_fwrite((void *) Y, out_integer, 1, fd) != 1) return DCD_BADWRITE;
  fio_write_int32(fd, out_integer);
  fio_write_int32(fd, out_integer);
  if (fio_fwrite((void *) Z, out_integer, 1, fd) != 1) return DCD_BADWRITE;
  fio_write_int32(fd, out_integer);

  /* update the DCD header information */
  fio_fseek(fd, NFILE_POS, FIO_SEEK_SET);
  fio_write_int32(fd, curframe);
  fio_fseek(fd, NSTEP_POS, FIO_SEEK_SET);
  fio_write_int32(fd, curstep);
  fio_fseek(fd, 0, FIO_SEEK_END);

  return DCD_SUCCESS;
}

/*
 * Write a header for a new dcd file
 * Input: fd - file struct opened for binary writing
 *        remarks - string to be put in the remarks section of the header.  
 *                  The string will be truncated to 70 characters.
 *        natoms, istart, nsavc, delta - see comments in read_dcdheader
 * Output: 0 on success, negative error code on failure.
 * Side effects: Header information is written to the dcd file.
 */
static int write_dcdheader(fio_fd fd, const char *remarks, int N, 
                    int ISTART, int NSAVC, double DELTA, int with_unitcell,
                    int charmm) {
  int out_integer;
  float out_float;
  char title_string[200];
  time_t cur_time;
  struct tm *tmbuf;
  char time_str[81];

  out_integer = 84;
  WRITE(fd, (char *) & out_integer, sizeof(int));
  strcpy(title_string, "CORD");
  WRITE(fd, title_string, 4);
  fio_write_int32(fd, 0);      /* Number of frames in file, none written yet   */
  fio_write_int32(fd, ISTART); /* Starting timestep                            */
  fio_write_int32(fd, NSAVC);  /* Timesteps between frames written to the file */
  fio_write_int32(fd, 0);      /* Number of timesteps in simulation            */
  fio_write_int32(fd, 0);      /* NAMD writes NSTEP or ISTART - NSAVC here?    */
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  if (charmm) {
    out_float = DELTA;
    WRITE(fd, (char *) &out_float, sizeof(float));
    if (with_unitcell) {
      fio_write_int32(fd, 1);
    } else {
      fio_write_int32(fd, 0);
    }
  } else {
    WRITE(fd, (char *) &DELTA, sizeof(double));
  }
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  if (charmm) {
    fio_write_int32(fd, 24); /* Pretend to be CHARMM version 24 */
  } else {
    fio_write_int32(fd, 0);
  }
  fio_write_int32(fd, 84);
  fio_write_int32(fd, 164);
  fio_write_int32(fd, 2);

  strncpy(title_string, remarks, 80);
  title_string[79] = '\0';
  WRITE(fd, title_string, 80);

  cur_time=time(NULL);
  tmbuf=localtime(&cur_time);
  strftime(time_str, 80, "REMARKS Created %d %B, %Y at %R", tmbuf);
  WRITE(fd, time_str, 80);

  fio_write_int32(fd, 164);
  fio_write_int32(fd, 4);
  fio_write_int32(fd, N);
  fio_write_int32(fd, 4);

  return DCD_SUCCESS;
}

static void *open_dcd_write(const char *path, const char *filetype, 
    int natoms) {
  dcdhandle *dcd;
  fio_fd fd;
  int rc;
  int istart, nsavc;
  double delta;
  int with_unitcell;
  int charmm;

  if (fio_open(path, FIO_WRITE, &fd) < 0) {
    printf("dcdplugin) Could not open file '%s' for writing\n", path);
    return NULL;
  }

  dcd = (dcdhandle *)malloc(sizeof(dcdhandle));
  memset(dcd, 0, sizeof(dcdhandle));
  dcd->fd = fd;

  istart = 0;             /* starting timestep of DCD file                  */
  nsavc = 1;              /* number of timesteps between written DCD frames */
  delta = 1.0;            /* length of a timestep                           */

  if (getenv("VMDDCDWRITEXPLORFORMAT") != NULL) {
    with_unitcell = 0;      /* no unit cell info */
    charmm = DCD_IS_XPLOR;  /* X-PLOR format */
    printf("dcdplugin) WARNING: Writing DCD file in X-PLOR format, \n");
    printf("dcdplugin) WARNING: unit cell information will be lost!\n");
  } else {
    with_unitcell = 1;      /* contains unit cell infor (Charmm format) */
    charmm = DCD_IS_CHARMM; /* charmm-formatted DCD file                */ 
    if (with_unitcell) 
      charmm |= DCD_HAS_EXTRA_BLOCK;
  }
 
  rc = write_dcdheader(dcd->fd, "Created by DCD plugin", natoms, 
                       istart, nsavc, delta, with_unitcell, charmm);

  if (rc < 0) {
    print_dcderror("write_dcdheader", rc);
    fio_fclose(dcd->fd);
    free(dcd);
    return NULL;
  }

  dcd->natoms = natoms;
  dcd->nsets = 0;
  dcd->istart = istart;
  dcd->nsavc = nsavc;
  dcd->with_unitcell = with_unitcell;
  dcd->charmm = charmm;
  dcd->x = (float *)malloc(natoms * sizeof(float));
  dcd->y = (float *)malloc(natoms * sizeof(float));
  dcd->z = (float *)malloc(natoms * sizeof(float));
  return dcd;
}


static int write_timestep(void *v, const molfile_timestep_t *ts) { 
  dcdhandle *dcd = (dcdhandle *)v;
  int i, rc, curstep;
  float *pos = ts->coords;
  double unitcell[6];
  unitcell[0] = unitcell[2] = unitcell[5] = 1.0f;
  unitcell[1] = unitcell[3] = unitcell[4] = 90.0f;

  /* copy atom coords into separate X/Y/Z arrays for writing */
  for (i=0; i<dcd->natoms; i++) {
    dcd->x[i] = *(pos++); 
    dcd->y[i] = *(pos++); 
    dcd->z[i] = *(pos++); 
  }
  dcd->nsets++;
  curstep = dcd->istart + dcd->nsets * dcd->nsavc;

  unitcell[0] = ts->A;
  unitcell[2] = ts->B;
  unitcell[5] = ts->C;
  unitcell[1] = sin((M_PI_2 / 90.0) * (90.0 - ts->gamma)); /* cosAB */
  unitcell[3] = sin((M_PI_2 / 90.0) * (90.0 - ts->beta));  /* cosAC */
  unitcell[4] = sin((M_PI_2 / 90.0) * (90.0 - ts->alpha)); /* cosBC */

  rc = write_dcdstep(dcd->fd, dcd->nsets, curstep, dcd->natoms, 
                     dcd->x, dcd->y, dcd->z,
                     dcd->with_unitcell ? unitcell : NULL,
                     dcd->charmm);
  if (rc < 0) {
    print_dcderror("write_dcdstep", rc);
    return MOLFILE_ERROR;
  }

  return MOLFILE_SUCCESS;
}

static void close_file_write(void *v) {
  dcdhandle *dcd = (dcdhandle *)v;
  fio_fclose(dcd->fd);
  free(dcd->x);
  free(dcd->y);
  free(dcd->z);
  free(dcd);
}

using namespace std;
molfile_timestep_t timestep;
void *v;
dcdhandle *dcd;
int natoms;

int Output_DCD_init(char*InputFileName){
  char DCD_FILE_NAME[256];
  int i;
  int last_dot_pos=-1;
  for(i=0;InputFileName[i]!=0;i++){
    if(InputFileName[i]=='.')last_dot_pos=i;
    DCD_FILE_NAME[i]=InputFileName[i];
  }
  if(last_dot_pos==-1)last_dot_pos=i;
  DCD_FILE_NAME[last_dot_pos]='.';
  DCD_FILE_NAME[last_dot_pos+1]='d';
  DCD_FILE_NAME[last_dot_pos+2]='c';
  DCD_FILE_NAME[last_dot_pos+3]='d';
  DCD_FILE_NAME[last_dot_pos+4]=0;

  natoms = 0;
  for(int type_id=0;type_id<Types.size();type_id++) {
    natoms+=Types[type_id].X.size();
  }
  v = open_dcd_write(DCD_FILE_NAME, "dcd", natoms);
  if (!v) {
    fprintf(stderr, "main) open_dcd_write failed for file %s\n", DCD_FILE_NAME);
    return 1;
  }
  dcd = (dcdhandle *)v;
  printf("main) file: %s\n", DCD_FILE_NAME);
  timestep.coords = (float *)malloc(3*sizeof(float)*natoms);
  timestep.A=Lx;
  timestep.B=Ly;
  timestep.C=Lz;
  timestep.alpha=90;
  timestep.beta =90;
  timestep.gamma=90;
  return 0;
}

int Output_DCD() {
  int index=-1;
  int type_id,bead_id;
  /*
  for(int type_id=0;type_id<Types.size();type_id++) {
    for(int bead_id=0;bead_id<Types[type_id].X.size();bead_id++){
      index++;timestep.coords[index]=Types[type_id].X[bead_id].x;
      index++;timestep.coords[index]=Types[type_id].X[bead_id].y;
      index++;timestep.coords[index]=Types[type_id].X[bead_id].z;
    }
  }*/
  for(int k=0;k<IO_Order_of_Beads.size();k++){
    type_id=IO_Order_of_Beads[k].x;
    bead_id=IO_Order_of_Beads[k].y;
    index++;timestep.coords[index]=Types[type_id].X[bead_id].x;
    index++;timestep.coords[index]=Types[type_id].X[bead_id].y;
    index++;timestep.coords[index]=Types[type_id].X[bead_id].z;
  }
  int res;
  res=write_timestep(dcd,&timestep);
  return res;
}

int Output_DCD_Close(){
  close_file_write(dcd);
  return 0;
}