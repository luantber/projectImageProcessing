#include "ift.h"

#define Imax 4095

/* Compute a gradient image from the multichannel image */

iftImage* MImageGradient(iftMImage* img, iftAdjRel* A) {
  iftImage* gradI = iftCreateImage(img->xsize, img->ysize, img->zsize);
  float Gmin = IFT_INFINITY_FLT, Gmax = IFT_INFINITY_FLT_NEG;
  iftFImage* fgradI = iftCreateFImage(img->xsize, img->ysize, img->zsize);

  for (int p = 0; p < img->n; p++) {
    iftVoxel u = iftMGetVoxelCoord(img, p);
    float value = 0.0;
    for (int i = 1; i < A->n; i++) {
      iftVoxel v = iftGetAdjacentVoxel(A, u, i);
      if (iftMValidVoxel(img, v)) {
        int q = iftMGetVoxelIndex(img, v);
        for (int b = 0; b < img->m; b++)
          value += (img->val[q][b] - img->val[p][b]) *
          (img->val[q][b] - img->val[p][b]);
      }
    }
    value = sqrtf(value) / A->n;
    if (value > Gmax) Gmax = value;
    if (value < Gmin) Gmin = value;
    fgradI->val[p] = value;
  }

  for (int p = 0; p < img->n; p++) {
    gradI->val[p] = (int)(Imax * (fgradI->val[p] - Gmin) / (Gmax - Gmin));
  }

  iftDestroyFImage(&fgradI);

  return(gradI);
}


iftFImage* MImageFGradient(iftMImage* img, iftAdjRel* A) {
  float Gmin = IFT_INFINITY_FLT, Gmax = IFT_INFINITY_FLT_NEG;
  iftFImage* fgradI = iftCreateFImage(img->xsize, img->ysize, img->zsize);

  for (int p = 0; p < img->n; p++) {
    iftVoxel u = iftMGetVoxelCoord(img, p);
    float value = 0.0;
    for (int i = 1; i < A->n; i++) {
      iftVoxel v = iftGetAdjacentVoxel(A, u, i);
      if (iftMValidVoxel(img, v)) {
        int q = iftMGetVoxelIndex(img, v);
        for (int b = 0; b < img->m; b++)
          value += (img->val[q][b] - img->val[p][b]) *
          (img->val[q][b] - img->val[p][b]);
      }
    }
    value = sqrtf(value) / A->n;
    if (value > Gmax) Gmax = value;
    if (value < Gmin) Gmin = value;
    fgradI->val[p] = value;
  }

  for (int p = 0; p < img->n; p++) {
    fgradI->val[p] = (Imax * (fgradI->val[p] - Gmin) / (Gmax - Gmin));
  }


  return(fgradI);
}

/* Compute a simple watershed transform on a gradient image */

iftImage* Watershed(iftMImage* img, iftLabeledSet* seeds) {
  iftImage* pathval = NULL, * label = NULL, * gradI = NULL;
  iftGQueue* Q = NULL;
  int         i, p, q, tmp;
  iftVoxel    u, v;
  iftLabeledSet* S = seeds;
  iftAdjRel* A = iftCircular(1.0);

  // Initialization

  pathval = iftCreateImage(img->xsize, img->ysize, img->zsize);
  label = iftCreateImage(img->xsize, img->ysize, img->zsize);
  gradI = MImageGradient(img, A);
  Q = iftCreateGQueue(Imax + 1, img->n, pathval->val);

  /* Initialize costs */

  while (S != NULL) {
    p = S->elem;
    label->val[p] = S->label;
    pathval->val[p] = 0;
    iftInsertGQueue(&Q, p);
    S = S->next;
  }

  for (p = 0; p < img->n; p++) {
    if (Q->L.elem[p].color == IFT_WHITE) { /* it is not seed */
      pathval->val[p] = IFT_INFINITY_INT;
      iftInsertGQueue(&Q, p);
    }
  }


  /* Propagate Optimum Paths by the Image Foresting Transform */

  while (!iftEmptyGQueue(Q)) {
    p = iftRemoveGQueue(Q);
    u = iftGetVoxelCoord(gradI, p);

    for (i = 1; i < A->n; i++) {
      v = iftGetAdjacentVoxel(A, u, i);

      if (iftValidVoxel(gradI, v)) {
        q = iftGetVoxelIndex(gradI, v);

        if (Q->L.elem[q].color != IFT_BLACK) {

          tmp = iftMax(pathval->val[p], gradI->val[q]);

          if (tmp < pathval->val[q]) {
            iftRemoveGQueueElem(Q, q);
            label->val[q] = label->val[p];
            pathval->val[q] = tmp;
            iftInsertGQueue(&Q, q);
          }
        }
      }
    }
  }

  iftDestroyAdjRel(&A);
  iftDestroyGQueue(&Q);
  iftDestroyImage(&gradI);
  iftDestroyImage(&pathval);

  return (label);
}

/* Compute a simple watershed transform on a float gradient image */

iftImage* FWatershed(iftMImage* img, iftLabeledSet* seeds) {
  iftFImage* pathval = NULL, * gradI = NULL;
  iftImage* label = NULL;
  iftFHeap* Q = NULL;
  int         i, p, q, tmp;
  iftVoxel    u, v;
  iftLabeledSet* S = seeds;
  iftAdjRel* A = iftCircular(1.0);

  // Initialization

  pathval = iftCreateFImage(img->xsize, img->ysize, img->zsize);
  label = iftCreateImage(img->xsize, img->ysize, img->zsize);
  gradI = MImageFGradient(img, A);
  Q = iftCreateFHeap(img->n, pathval->val);

  /* Initialize costs */

  while (S != NULL) {
    p = S->elem;
    label->val[p] = S->label;
    pathval->val[p] = 0;
    iftInsertFHeap(Q, p);
    S = S->next;
  }

  for (p = 0; p < img->n; p++) {
    if (Q->color[p] == IFT_WHITE) { /* it is not seed */
      pathval->val[p] = IFT_INFINITY_FLT;
      iftInsertFHeap(Q, p);
    }
  }

  /* Propagate Optimum Paths by the Image Foresting Transform */

  while (!iftEmptyFHeap(Q)) {
    p = iftRemoveFHeap(Q);
    u = iftFGetVoxelCoord(gradI, p);

    for (i = 1; i < A->n; i++) {
      v = iftGetAdjacentVoxel(A, u, i);

      if (iftFValidVoxel(gradI, v)) {
        q = iftFGetVoxelIndex(gradI, v);

        if (Q->color[q] != IFT_BLACK) {

          tmp = iftMax(pathval->val[p], gradI->val[q]);

          if (tmp < pathval->val[q]) {
            label->val[q] = label->val[p];
            pathval->val[q] = tmp;
            iftGoUpFHeap(Q, Q->pos[q]);
          }
        }
      }
    }
  }

  for (int i = 0; i < label->n; i++) {
    label->val[i] = label->val[i] * 10;
  }


  iftDestroyAdjRel(&A);
  iftDestroyFHeap(&Q);
  iftDestroyFImage(&gradI);
  iftDestroyFImage(&pathval);

  return (label);
}



/* Extracts the basename of an image file */

char* Basename(char* path) {
  char* basename = iftBasename(path);
  iftSList* slist = iftSplitString(basename, "/");
  strcpy(basename, slist->tail->elem);
  iftDestroySList(&slist);
  return(basename);
}


int main(int argc, char* argv[]) {
  timer* tstart;
  char       filename[200];

  if (argc != 4)
    iftError("Usage: watershed <...>\n"
      "[1] input image .png \n"
      "[2] labeled seed set .txt \n"
      "[3] output folder\n",
      "main");

  tstart = iftTic();

  /* Read the input image and convert it to the Lab Color Space */

  iftImage* img = iftReadImageByExt(argv[1]);
  char* basename = Basename(argv[1]);
  iftMImage* mimg = NULL;

  mimg = iftImageToMImage(img, LABNorm2_CSPACE);

  /* Create output folder and read labeled seeds */

  iftMakeDir(argv[3]);
  iftLabeledSet* S = iftReadSeeds(img, argv[2]);

  /* Compute watershed segmentation */

  iftImage* label = FWatershed(mimg, S);

  /* Draw object borders and save the output images */

  iftAdjRel* A = iftCircular(1.0);
  iftColor RGB, YCbCr;
  RGB.val[0] = 255;
  RGB.val[1] = 0;
  RGB.val[2] = 255;
  YCbCr = iftRGBtoYCbCr(RGB, 255);
  iftDrawBorders(img, label, A, YCbCr, A);
  iftDestroyAdjRel(&A);

  sprintf(filename, "%s/%s_label.png", argv[3], basename);
  iftWriteImageByExt(label, filename);


  sprintf(filename, "%s/%s_segm.png", argv[3], basename);
  iftWriteImageByExt(img, filename);


  iftDestroyImage(&img);
  iftDestroyImage(&label);
  iftDestroyMImage(&mimg);
  iftFree(basename);
  iftDestroyLabeledSet(&S);

  printf("Done ... %s\n", iftFormattedTime(iftCompTime(tstart, iftToc())));
  return (0);
}
