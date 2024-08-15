#ifndef SRC_S21_MATRIX_H_
#define SRC_S21_MATRIX_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define OK 0
#define INCORRECT_MATRIX 1
#define CALCULATION_ERROR 2
#define CALC_ERROR 2
#define ALLOC_FAIL 3

#define SUCCESS 1
#define FAILURE 0

// 0 - OK
// 1 - Ошибка, некорректная матрица
// 2 - Ошибка вычисления (несовпадающие размеры матриц; матрица, для которой
// нельзя провести вычисления и т.д.)

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

int s21_create_matrix(int rows, int columns, matrix_t *result);
int s21_eq_matrix(matrix_t *A, matrix_t *B);
void s21_remove_matrix(matrix_t *A);
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int check_error(matrix_t *matrix);
int s21_transpose(matrix_t *A, matrix_t *result);
int s21_calc_complements(matrix_t *A, matrix_t *result);
int s21_matrix_minor(matrix_t *A, matrix_t *result, int i, int j);
int s21_determinant(matrix_t *A, double *result);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);
double s21_get_determinant(matrix_t *A);
void s21_get_matrix(int row, int col, matrix_t *A, matrix_t *result);

#endif