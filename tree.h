#include "stdlib.h"
typedef struct Nodo {

  struct Nodo* next;
  float* value;

} Nodo;


typedef struct Tree {
  Nodo** roots;
  Nodo** lasts;

  float** acumular;
  int* acumular_size;

  int n;
  int bandas;
} Tree;


Tree* initTree(int n, int bandas) {
  Tree* tree = (Tree*)malloc(sizeof(Tree));
  tree->roots = (Nodo**)malloc(sizeof(Nodo*) * n);
  tree->lasts = (Nodo**)malloc(sizeof(Nodo*) * n);

  tree->acumular = (float**)malloc(sizeof(float*) * n);
  for (int i = 0; i < n; i++) {
    tree->acumular[i] = (float*)malloc(sizeof(float) * bandas);
  }

  tree->acumular_size = (int*)malloc(sizeof(int) * n);

  for (int i = 0; i < n; i++) {
    tree->roots[i] = NULL;
    tree->lasts[i] = NULL;

    for (int j = 0; j < bandas; j++) {
      tree->acumular[i][j] = 0;
    }

    tree->acumular_size[i] = 0;

  }
  tree->n = n;
  tree->bandas = bandas;
  return tree;
}

void appendTree(Tree* tree, int r, float* value) {
  Nodo* newNode = (Nodo*)malloc(sizeof(Nodo));
  newNode->value = value;
  newNode->next = NULL;

  if (tree->roots[r] == NULL) {
    tree->roots[r] = newNode;
    tree->lasts[r] = newNode;
  }
  else {
    tree->lasts[r]->next = newNode;
    tree->lasts[r] = newNode;
  }

  for (int i = 0; i < tree->bandas; i++) {

    tree->acumular[r][i] += value[i];

  }
  tree->acumular_size[r] += 1;


}

float* averageTree(Tree* tree, int r) {
  // printf("averaging tree\n");


  float* sum = (float*)malloc(sizeof(float) * tree->bandas);

  for (int i = 0; i < tree->bandas; i++) {

    sum[i] = tree->acumular[r][i] / tree->acumular_size[r];
  }

  return sum;

}

float normaRestarArrays(float* a, float* b, int m) {
  // printf("restando tree\n");
  float sum = 0;
  for (int i = 0; i < m; i++) {
    sum += fabs(a[i] - b[i]);
  }
  return sum;
}

void freeTree(Tree* tree) {
  for (int i = 0; i < tree->n; i++) {
    Nodo* current = tree->roots[i];
    while (current != NULL) {
      Nodo* next = current->next;
      free(current);
      current = next;
    }
  }
  free(tree->roots);
  free(tree->lasts);
  for (int i = 0; i < tree->n; i++) {
    free(tree->acumular[i]);
  }
  free(tree->acumular);
  free(tree->acumular_size);
  free(tree);
}