#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include "tinyexpr.h"
#define A_EPS 1e-12
#define H 1e-12



double solve_f(int arg_num, char *problem, int dx, double xForDx, double arr_x[])
{

    double c = arr_x[dx];
    arr_x[dx] = xForDx;

    const char *expression = problem;

    te_variable* vars = malloc((arg_num - 2) * sizeof(te_variable));
    char* ch[arg_num - 2][4];
    for (int i = 0; i < arg_num - 2; i++)
    {
        sprintf(ch[i], "%c%d", 'x', i + 1);
    }
    for (int i = 0; i < arg_num - 2; i++)
    {
        vars[i].name = ch[i];
        vars[i].address = &arr_x[i];
    }

    te_expr *n = te_compile(expression, vars, arg_num - 2, 0);
    if (n)
        {
            const double rf = te_eval(n);
            te_free(n);
            arr_x[dx] = c;
            te_free(vars);
            return rf;
        }

}


double diff_dfdx(int arg_num, char *problem, int dx, double arr_x[])
{
    return (solve_f(arg_num, problem, dx, arr_x[dx] + H, arr_x) - solve_f(arg_num, problem, dx, arr_x[dx] - H, arr_x)) / (2.0 * H);
}

double solve_dfdx(int arg_num, char *problem, int dx, double arr_x[])
{
    return diff_dfdx(arg_num, problem, dx, arr_x);
}


void solve_main(int size, char *problem[][50], double arr_solve_f[], double arr_x[])
{

    double arr_m[size][size];

    //Creating identity matrix:
    double arr_idm[size][size];
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            if (i == j)
                arr_idm[i][j] = 1;
            else
                arr_idm[i][j] = 0;
        }
    }

    //Creating common matrix
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            arr_m[i][j] = solve_dfdx(size + 2, problem[i], j, arr_x);
        }
    }


    //Evaluating inverse matrix:
    double r, res;
    for (int k = 0; k < size; k++)
    {
        r = 1/arr_m[k][k];
        for (int j = 0; j < size; j++)
        {
            arr_m[k][j] *= r;
            arr_idm[k][j] *= r;
        }
        for (int i = k + 1; i < size; i++)
        {
            res = arr_m[i][k];
            for (int z = 0; z < size; z++)
            {
                arr_m[i][z] -= arr_m[k][z] * res;
                arr_idm[i][z] -= arr_idm[k][z] * res;
            }
        }
    }
    for (int k = size - 1; k >= 0; k--)
    {
        for (int i = k - 1; i >= 0; i--)
        {
            res = arr_m[i][k];
            for (int z = size - 1; z >= 0; z--)
            {
                arr_m[i][z] = arr_m[i][z] - arr_m[k][z] * res;
                arr_idm[i][z] = arr_idm[i][z] - arr_idm[k][z] * res;
            }
        }
    }


    for (int i = 0; i < size; i++)
    {
        arr_solve_f[i] = 0;
        for (int j = 0; j < size; j++)
        {

            arr_solve_f[i] += arr_idm[i][j] * solve_f(size + 2, problem[j], 0, arr_x[0], arr_x);

        }
        arr_solve_f[i] = arr_x[i] - arr_solve_f[i];
    }

}






int main(void)
{
    int exp_num;
    scanf("%d", &exp_num);

    double arr_x[exp_num]; 
    for (int i = 0; i < exp_num; i++)
    {
        double x;
        scanf("%lf", &x);
        arr_x[i] = x;
    }

    for (int i = 0; i < exp_num; i++)
    {
        printf("x[%d] = %lf ", i, arr_x[i]);
    }

    char *problem[exp_num][50];
    for (int i = 0; i < exp_num; i++)
    {
        scanf("%s", &problem[i]);
    }


    double arr_solve_f[exp_num];


    if (exp_num == 1)
    {
        solve_f(exp_num + 2, problem[0], 0, arr_x[0], arr_x);
        solve_dfdx(exp_num + 2, problem[0], 0, arr_x);
        while (fabs(solve_f(exp_num + 2, problem[0], 0, arr_x[0], arr_x)) > A_EPS)
        {
            arr_x[0] -= solve_f(exp_num + 2, problem[0], 0, arr_x[0], arr_x) / solve_dfdx(exp_num + 2, problem[0], 0, arr_x);
        }
        printf("x1 = %lf\n", arr_x[0]);
        return 0;
    }
    else
    {
        solve_f(exp_num + 2, problem[0], 0, arr_x[0], arr_x);
        solve_dfdx(exp_num + 2, problem[1], 0, arr_x);
        while (fabs(solve_f(exp_num + 2, problem[1], 0, arr_x[0], arr_x)) > A_EPS)
        {
            solve_main(exp_num, problem, arr_solve_f, arr_x);
            for (int i = 0; i < exp_num; i++)
            {
                arr_x[i] = arr_solve_f[i];
            }
        }

        for (int i = 0; i < exp_num; i++)
            printf("x%d = %lf ", i + 1, arr_x[i]);

    }


    return 0;
}
