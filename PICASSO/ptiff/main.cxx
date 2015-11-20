/* PICASSO TIFF viewer by Ethan S. Miller (esmiller@uiuc.edu) */

#include <FL/Fl_PNG_Image.H>
#include <FL/Fl_JPEG_Image.H>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>

#include "ptiffUI.h"
#include "Fl_TIFF_Image.H"

using namespace std;

int main (int argc, char **argv) {
  ptiffUI *UI = new ptiffUI;
  UI->window_->show();

  fl_register_images();
  Fl_Shared_Image::add_handler(Fl_TIFF_Image::check);

  return Fl::run();
}
