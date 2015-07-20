/*
** Include file to be used in conjuntion with the volume set (Volumes)
** data type.
**
** Author:
**  G.Lohmann, <lohmann@cns.mpg.de>, Jan. 1996
*/

#include <stdio.h>

#define MAXHASHLEN 1024  /* max length of hash table */

typedef struct VTrackStruct {
  short band;
  short row; 
  short col;
  short length;
  struct VTrackStruct *next;
  struct VTrackStruct *previous;
} *VTrack, VTrackRec;


typedef struct VBucketStruct {
  short  ntracks;      /* number of tracks in one hashtable bucket */
  VTrack first;        /* ptr to first track in bucket             */
  VTrack last;         /* ptr to last track in bucket              */
} *VBucket, VBucketRec;

typedef struct VolumeStruct {
  short label;
  short nbands;
  short nrows;
  short ncolumns;
  short nbuckets;      /* length of hash table (number of buckets) */
  int   ntracks;       /* total number of tracks in all buckets   */
  VBucket bucket;      /* ptrs to buckets      */
  struct VolumeStruct *next;
} VolumeRec, *Volume;


typedef struct V_VolumesRec {
  VAttrList attributes;
  short nvolumes;        /* number of volumes in list       */
  short nbands;
  short nrows;
  short ncolumns;
  Volume first;          /* ptr to first volume in list     */
} VolumesRec;

#define VolumesAttrList(volumes) ((volumes)->attributes)
#define VolumesNum(volumes) ((volumes)->nvolumes)
#define VolumesNBands(volumes) ((volumes)->nbands)
#define VolumesNRows(volumes) ((volumes)->nrows)
#define VolumesNColumns(volumes) ((volumes)->ncolumns)
#define VolumesNVolumes(volumes) ((volumes)->ntracks)
/*
#define VolumesNTracks(volumes) ((volumes)->ntracks)
*/

#define VolumeNBands(volume) ((volume)->nbands)
#define VolumeNRows(volume) ((volume)->nrows)
#define VolumeNColumns(volume) ((volume)->ncolumns)
#define VolumeNBuckets(volume) ((volume)->nbuckets)
#define VolumeNTracks(volume) ((volume)->ntracks)
#define VolumeLabel(volume) ((volume)->label)
#define VFirstVolume(volumes) ((volumes)->first)
#define VNextVolume(volume) ((volume)->next)
#define VolumeExists(volume) ((volume) != NULL)

#define VTrackLength(track) ((track)->length)
#define VTrackExists(track) ((track) != NULL)
#define VFirstTrack(volume,i) ((volume)->bucket[(i)].first)
#define VNextTrack(track) ((track)->next)
#define VPreviousTrack(track) ((track)->previous)

#define VolumesAttr     "volumes"
#define VolNVolumesAttr "nvolumes"
#define VolNTracksAttr  "ntracks"
#define VolNBandsAttr   "nbands"
#define VolNRowsAttr    "nrows"
#define VolNColumnsAttr "ncolumns"

extern Volumes VCreateVolumes(short,short,short);
extern Volumes VCopyVolumes(Volumes);
extern void VDestroyVolumes(Volumes);
extern VBoolean VWriteVolumes(FILE *, VAttrList, int, Volumes *);

extern int VReadVolumes(FILE *, VAttrList *, Volumes **);

extern Volume VCreateVolume(short,short,short,short,short);
extern Volume VCopyVolume(Volume);
extern void VAddVolume(Volumes, Volume);
extern void AddTrack(Volume,VTrack);

extern double VolumeBorderSize(Volume);
extern VBoolean VolumeBorder(Volume,short,short,short);
extern VTrack VolumeGetTrack(Volume,short,short,short);
extern VBoolean VolumeInside(Volume,short,short,short);
extern double VolumeRadius(Volume,double *);

/*
** hash function
*/
#define VolumeHash(nbands, b, r, len) (((b) * (nbands) + (r)) % (len))
