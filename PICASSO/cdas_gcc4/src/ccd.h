/*
 * For C++ compilers, use extern "C"
 */

#ifdef __cplusplus
extern "C" {
#endif



int tcl_read_image(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);
int tcl_write_image(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);
int tcl_write_cimage(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);
int tcl_write_fimage(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);
int tcl_write_dimage(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);
int tcl_write_zimage(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);
int tcl_write_simage(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);
int tcl_shmem_image(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);
int tcl_show_image(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);
int tcl_write_calibrated(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);
int tcl_store_calib(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);
int tcl_list_buffers(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);
int tcl_set_biascols(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);
int tcl_set_biasrows(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);
int writeimage(void *src_buffer, char *filename, short nx, short by, short ibx,short iby);
int readimage(char *buffer_name, char *filename);
int calibrateimage(void *src_buffer, char *filename);
int printerror(int status);
int biascorrimage(void *src_buffer, char *filename, short nx, short ny);
int biassubtract(unsigned short *src,float *dest, int nx, int ny);
int ibiassubtract(unsigned short *src,unsigned short *dest, int nx, int ny);
int write_buffered_image(void *buffer_name, char *filename);
int calculate_flatfield(float *image,float *work);
int calculate_dark(float *image);
int calculate_zero(float *image);
void determine_fbiascols(int buffernum);
void determine_ibiascols(int buffernum);
int calculate_skyflat(float *image, float *work);
float calculate_median(float *data, int n);
int image_i2tof(unsigned short *src,float *dest,int n);
int image_rawi2tof(unsigned short *src,float *dest,int n);
void subtractzero(float *rimg,float *temp,int nelements);
void subtractdark(float *rimg,float *temp,int nelements);
void divide(float *rimg,float *temp,int nelements);
int shmmapimage(void *buffer_name);
int disp_init (int w, int h, int fbconfig, int frame);
void converttobyte(float *src, unsigned char *dest,int n);

typedef void *PDATA;

PDATA CCD_locate_buffer(char *name, int idepth, short imgcols, short imgrows, short hbin, short vbin);

int   CCD_free_buffer();
int   CCD_locate_buffernum(char *name);


/*
 * end block for C++
 */

#ifdef __cplusplus
}
#endif

