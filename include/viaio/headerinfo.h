typedef struct V_ImageInfo {
  long nbands;                       /* number of bands */
  long nrows;                        /* number of rows */
  long ncolumns;                     /* number of columns */
  VRepnKind repn;                    /* representation of pixel values */
  long pixelSize;                     /* number of bytes per pixel */
  size_t data;                       /* address of first pixel */
  size_t length;                     /* number of bytes in image */
  VShort ori_nbands;                 /* info needed for compressed format */
  VShort ori_nrows;
  VShort ori_ncolumns;
  VShort top_margin;
  VShort left_margin;
  VString patient;
  VString modality;
  VString angle;
  VString voxel;
  VString name;
  VString fixpoint;
  VString ca;
  VString cp;
  VString location;
  VString extent;
  VString orientation;
  VString talairach;
  VString MPIL_vista_0;
  size_t offsetHdr;                    /* length of ASCII header in bytes */
  long   spPH;
  long   spPG;
  long   subjects;
  long   ntimesteps;
  VDouble df;
  VDouble norm_mean;
  VDouble norm_sig;
  long repetition_time;
} VImageInfo;

