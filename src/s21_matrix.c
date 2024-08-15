#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int res = OK;
  if (rows <= 0 || columns <= 0 || result == NULL) {
    return INCORRECT_MATRIX;
  }
  result->rows = rows;
  result->columns = columns;
  result->matrix = (double **)calloc(rows, sizeof(double *));
  if (result->matrix == NULL) {
    return INCORRECT_MATRIX;
  }
  for (int i = 0; i < rows; i++) {
    result->matrix[i] = (double *)calloc(columns, sizeof(double));
    if (result->matrix[i] == NULL) {
      res = INCORRECT_MATRIX;
      for (int j = 0; j < i; j++) {
        free(result->matrix[j]);  // освобождаем уже выделенную память
      }
      free(result->matrix);
      result->matrix = NULL;
      break;
    }
    memset(result->matrix[i], 0,
           columns * sizeof(double));  // инициализируем значениями 0
  }
  if (res != OK) {
    result->rows = 0;
    result->columns = 0;
  }
  return res;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int res = SUCCESS;
  if (check_error(A) != OK || check_error(B) != OK ||
      A->columns != B->columns || A->rows != B->rows) {
    res = FAILURE;
  } else {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) >= 1e-7) res = FAILURE;
      }
    }
  }

  return res;
}

void s21_remove_matrix(matrix_t *A) {
  if (A->matrix != NULL) {
    for (int i = 0; i < (A->rows); i++) {
      free(A->matrix[i]);
      A->matrix[i] = NULL;
    }
    free(A->matrix);
    A->matrix = NULL;
  }
  A->columns = 0;
  A->rows = 0;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = OK;
  if (A != NULL && B != NULL) {
    if (A->rows != B->rows || A->columns != B->columns) {
      res = CALCULATION_ERROR;
    } else {
      res = s21_create_matrix(A->rows, A->columns, result);
      if (res == 0) {
        for (int i = 0; i < A->rows; i++) {
          for (int j = 0; j < A->columns; j++) {
            result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
          }
        }
      }
    }

  } else {
    res = INCORRECT_MATRIX;
  }
  return res;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = OK;
  if (A != NULL && B != NULL) {
    if (A->rows != B->rows || A->columns != B->columns) {
      res = CALCULATION_ERROR;
    } else {
      res = s21_create_matrix(A->rows, A->columns, result);
      if (res == 0) {
        for (int i = 0; i < A->rows; i++) {
          for (int j = 0; j < A->columns; j++) {
            result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
          }
        }
      }
    }

  } else {
    res = INCORRECT_MATRIX;
  }
  return res;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int res = OK;
  int error = check_error(A);
  if (error == 0) {
    res = s21_create_matrix(A->rows, A->columns, result);
    if (res == 0) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] * number;
        }
      }
    }
  } else {
    res = INCORRECT_MATRIX;
  }
  return res;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = OK;
  int error_a = check_error(A);
  int error_b = check_error(B);

  if (error_a == OK && error_b == OK) {
    if (A->columns == B->rows) {
      res = s21_create_matrix(A->rows, B->columns, result);
      // явная инициализация нулями
      for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->columns; j++) {
          result->matrix[i][j] = 0.0;
        }
      }
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < B->columns; j++) {
          for (int k = 0; k < B->rows; k++) {
            result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
          }
        }
      }
    } else {
      res = CALCULATION_ERROR;
    }
  } else {
    res = INCORRECT_MATRIX;
  }
  return res;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int res = OK;
  if (check_error(A) == 0) {
    res = s21_create_matrix(A->columns, A->rows, result);
    if (res == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int k = 0; k < A->columns; k++) {
          result->matrix[k][i] = A->matrix[i][k];
        }
      }
    }
  } else {
    res = INCORRECT_MATRIX;
  }
  return res;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int res = OK;
  if (check_error(A) == 0) {
    if (A->rows == A->columns) {
      s21_create_matrix(A->columns, A->rows, result);
      matrix_t minor = {0};
      s21_create_matrix(A->columns - 1, A->rows - 1, &minor);
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          s21_get_matrix(i, j, A, &minor);
          double complement = 0.0;
          s21_determinant(&minor, &complement);
          result->matrix[i][j] = pow(-1, i + j) * complement;
        }
      }
      s21_remove_matrix(&minor);
    } else {
      res = CALCULATION_ERROR;
    }
  } else {
    res = INCORRECT_MATRIX;
  }
  return res;
}

int s21_determinant(matrix_t *A, double *result) {
  *result = 0.0;
  int res = OK;
  if (check_error(A) == OK) {
    if (A->rows == A->columns) {
      *result = s21_get_determinant(A);
    } else {
      res = CALCULATION_ERROR;
    }
  } else {
    res = INCORRECT_MATRIX;
  }
  return res;
}

double s21_get_determinant(matrix_t *A) {
  double res = 0.0;
  if (A->rows == 1) {
    res = A->matrix[0][0];
  } else {
    matrix_t tmp = {0};
    s21_create_matrix(A->rows - 1, A->columns - 1, &tmp);
    for (int i = 0; i < A->columns; i++) {
      s21_get_matrix(0, i, A, &tmp);
      if (i % 2) {
        res -= A->matrix[0][i] * s21_get_determinant(&tmp);
      } else {
        res += A->matrix[0][i] * s21_get_determinant(&tmp);
      }
    }
    s21_remove_matrix(&tmp);
  }
  return res;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int res = INCORRECT_MATRIX;
  if (check_error(A) == OK) {
    if (A->rows == A->columns) {  // Добавляем проверку на квадратность матрицы
      double det = 0.0;
      s21_determinant(A, &det);
      if (det != 0) {
        matrix_t tmp_calc = {0};
        res = s21_calc_complements(A, &tmp_calc);
        if (res == OK) {
          matrix_t tmp_trans = {0};
          res = s21_transpose(&tmp_calc, &tmp_trans);
          if (res == OK) {
            s21_mult_number(&tmp_trans, 1 / det, result);
          }
          s21_remove_matrix(&tmp_trans);
        }
        s21_remove_matrix(&tmp_calc);
      } else {
        res = CALCULATION_ERROR;
      }
    } else {
      res = CALCULATION_ERROR;  // Возвращаем CALCULATION_ERROR, если матрица не
                                // квадратная
    }
  }
  return res;
}

int check_error(matrix_t *matrix) {
  int res = OK;
  if (matrix == NULL || matrix->matrix == NULL || matrix->rows <= 0 ||
      matrix->columns <= 0) {
    res = INCORRECT_MATRIX;
  } else {
    res = OK;
  }
  return res;
}

// выполняет операцию удаления строки
// и столбца с заданными индексами из матрицы A и
// записывает результат в результирующую матрицу minor.
void s21_get_matrix(int row, int col, matrix_t *A, matrix_t *minor) {
  int minor_row = 0;
  int minor_col = 0;

  for (int i = 0; i < A->rows; i++) {
    if (i == row) {
      continue;
    }
    minor_col = 0;
    for (int j = 0; j < A->columns; j++) {
      if (j == col) {
        continue;
      }
      minor->matrix[minor_row][minor_col] = A->matrix[i][j];
      minor_col++;
    }
    minor_row++;
  }
}
