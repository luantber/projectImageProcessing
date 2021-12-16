#include "ift.h"
#include "tree.h"
#define Imax 4095

/* Compute a gradient image from the multichannel image */
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

/* Compute a simple watershed transform on a float gradient image */
iftImage* FWatershed(iftMImage* img, iftLabeledSet* S) {

  // Cost function, Label, and 
  iftFImage* C = iftCreateFImage(img->xsize, img->ysize, img->zsize);
  iftImage* L = iftCreateImage(img->xsize, img->ysize, img->zsize);

  // Heap
  iftFHeap* Q = iftCreateFHeap(img->n, C->val);

  //Adjancecy relation
  iftAdjRel* A = iftCircular(1.0);

  // Gradient
  iftFImage* gradI = MImageFGradient(img, A);

  // Voxel coordinates
  iftVoxel    u_p, v_p;

  // creates a R(p) = q 
  int* R = (int*)malloc(img->n * sizeof(int));

  // creates a P(p) = q  // predecesor
  int* P = (int*)malloc(img->n * sizeof(int));

  int n_labels = iftNumberOfLabels(S);

  Tree* tree = initTree(img->n, img->m);

  // Initialize cost 
  while (S != NULL) {
    int p = S->elem;
    L->val[p] = S->label;
    C->val[p] = 0;

    // Dynamic
    R[p] = p;
    P[p] = -1;

    iftInsertFHeap(Q, p);
    S = S->next;
  }

  for (int p = 0; p < img->n; p++) {
    if (Q->color[p] == IFT_WHITE) { /* it is not seed */
      C->val[p] = IFT_INFINITY_FLT;

      //Dynamic
      R[p] = p;
      P[p] = -1;

      iftInsertFHeap(Q, p);
    }
  }


  /* Propagate Optimum Paths by the Image Foresting Transform */
  while (!iftEmptyFHeap(Q)) {

    int p = iftRemoveFHeap(Q);
    u_p = iftFGetVoxelCoord(C, p);

    appendTree(tree, R[p], img->val[p]);

    for (int i = 1; i < A->n; i++) {
      v_p = iftGetAdjacentVoxel(A, u_p, i);

      if (iftFValidVoxel(gradI, v_p)) {
        int q = iftFGetVoxelIndex(gradI, v_p);

        if (Q->color[q] != IFT_BLACK) {

          float w0 = gradI->val[q];

          // Dynamic by Average
          float* mu = averageTree(tree, R[p]);
          float w1 = normaRestarArrays(mu, img->val[q], img->m);


          float tmp = iftMax(C->val[p], w1);
          free(mu);

          if (tmp < C->val[q]) {

            C->val[q] = tmp;
            L->val[q] = L->val[p];

            // Dynamic
            R[q] = R[p];
            P[q] = p;

            iftGoUpFHeap(Q, Q->pos[q]);

          }

        }
      }
    }
  }

  for (int i = 0; i < L->n; i++) {
    L->val[i] = L->val[i] * 25;
  }

  iftDestroyAdjRel(&A);
  iftDestroyFHeap(&Q);
  iftDestroyFImage(&gradI);
  iftDestroyFImage(&C);

  free(R);
  free(P);
  freeTree(tree);

  return L;
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
