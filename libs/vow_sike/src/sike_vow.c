#include <stdio.h>
#include <time.h>
#include <string.h>
#include "bintree.h"
#include "sidh_vow_base.h"
#include "vow.h"
#include "memory.h"
#include "state.h"
#include "sike_vow.h"
#define BRANCHING_FACTOR 4

uint64_t *starting_index_for_each_level;
uint64_t max_index_w;

shared_container right_precomputation_tree;
shared_container E_left;
shared_container E_right;

void GetTwoIsogenyWithXneZ(const point_proj_t R, point_proj_t A24)
{
    fp2sqr_mont(R->X, A24->X);
    fp2sqr_mont(R->Z, A24->Z);
    fp2sub(A24->Z, A24->X, A24->X);
}

void EvalTwoIsogenyWithXneZ(const point_proj_t R, point_proj_t P)
{
    f2elm_t T0, T1, T2;

    fp2add(P->X, P->Z, T0);
    fp2sub(R->Z, R->X, T1);
    fp2mul_mont(T0, T1, T0);
    fp2sub(P->X, P->Z, T1);
    fp2add(R->X, R->Z, T2);
    fp2mul_mont(T1, T2, T1);
    fp2sub(T1, T0, T2);
    fp2add(T1, T0, T1);
    fp2mul_mont(P->X, T2, P->X);
    fp2mul_mont(P->Z, T1, P->Z);
}

static void print_curve_parameter(CurveAndPointsSIDH curve)
{
    f2elm_t a24;
    f2elm_t xp_standard;
    f2elm_t xq_standard;
    f2elm_t xpq_standard;
    from_mont(curve.a24[0], a24[0]);
    from_mont(curve.a24[1], a24[1]);
    printf("a24:\n");
    print_f2elm(a24);

    from_mont(curve.xp[0], xp_standard[0]);
    from_mont(curve.xp[1], xp_standard[1]);
    printf("xp:\n");
    print_f2elm(xp_standard);

    from_mont(curve.xq[0], xq_standard[0]);
    from_mont(curve.xq[1], xq_standard[1]);
    printf("xq:\n");
    print_f2elm(xq_standard);

    from_mont(curve.xpq[0], xpq_standard[0]);
    from_mont(curve.xpq[1], xpq_standard[1]);
    printf("xpq:\n");
    print_f2elm(xpq_standard);
}

void LadderThreePtSIKE(point_proj_t R, const f2elm_t a24, const f2elm_t xp, const f2elm_t xq, const f2elm_t xpq, const unsigned char *m, unsigned long nbits_state)
{ // Non-constant time version that depends on size of k
    point_proj_t R0 = {0}, R2 = {0};
    unsigned int i, bit;

    fp2copy(xp, R->X);
    fp2copy(xq, R0->X);
    fp2copy(xpq, R2->X);
    fpcopy((digit_t *)&Montgomery_one, (digit_t *)R0->Z);
    fpzero((digit_t *)(R0->Z)[1]);
    fpcopy((digit_t *)&Montgomery_one, (digit_t *)R->Z);
    fpzero((digit_t *)(R->Z)[1]);
    fpcopy((digit_t *)&Montgomery_one, (digit_t *)R2->Z);
    fpzero((digit_t *)(R2->Z)[1]);

    for (i = 0; i < nbits_state; i++)
    { // Ignore 3 lsb's c,b
        bit = (m[i >> 3] >> (i & 0x07)) & 1;

        if (bit)
        {
            xDBLADD_SIDH(R0, R, R2->X, a24);
            fp2mul_mont(R->X, R2->Z, R->X);
        }
        else
        {
            xDBLADD_SIDH(R0, R2, R->X, a24);
            fp2mul_mont(R2->X, R->Z, R2->X);
        }
    }
}

void IsogenyWithPoints(point_proj_t A24, point_proj_t XP, point_proj_t XQ, point_proj_t XPQ,
                       const f2elm_t a24, const point_proj_t kernel_point, const point_proj_t XQinp,
                       const unsigned long *strat, const unsigned long lenstrat)
{
    point_proj_t R, pts[MAX_INT_POINTS_ALICE];
    f2elm_t coeff[3];
    unsigned long i, row, index = 0, ii = 0, m, npts = 0, pts_index[MAX_INT_POINTS_ALICE];

    fp2copy(kernel_point->X, R->X);
    fp2copy(kernel_point->Z, R->Z);

    fp2copy(XQinp->X, XQ->X);
    fp2copy(XQinp->Z, XQ->Z);

    fp2copy(a24, A24->X);
    fpcopy((digit_t *)&Montgomery_one, (digit_t *)A24->Z);
    fpzero((digit_t *)(A24->Z)[1]);

    /* All steps except the first */
    for (row = 1; row < lenstrat + 1; row++)
    {
        while (index < lenstrat + 1 - row)
        {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            pts_index[npts++] = index;
            m = strat[ii++];
            xDBLe(R, R, A24->X, A24->Z, (int)(2 * m));
            index += m;
        }
        GetFourIsogenyWithKernelXneZ(R, A24->X, A24->Z, coeff);
        EvalFourIsogenyWithKernelXneZ(XP, coeff);
        EvalFourIsogenyWithKernelXneZ(XQ, coeff);
        EvalFourIsogenyWithKernelXneZ(XPQ, coeff);

        for (i = 0; i < npts; i++)
        {
            EvalFourIsogenyWithKernelXneZ(pts[i], coeff);
        }

        fp2copy(pts[npts - 1]->X, R->X);
        fp2copy(pts[npts - 1]->Z, R->Z);
        index = pts_index[npts - 1];
        npts -= 1;
    }

    GetFourIsogenyWithKernelXneZ(R, A24->X, A24->Z, coeff);
    EvalFourIsogenyWithKernelXneZ(XP, coeff);
    EvalFourIsogenyWithKernelXneZ(XQ, coeff);
    EvalFourIsogenyWithKernelXneZ(XPQ, coeff);
}

static void copy_CurveAndPointsSIDH(CurveAndPointsSIDH *dst, CurveAndPointsSIDH *src)
{
    // fp2copy(src, dst)
    fp2copy(src->a24, dst->a24);
    fp2copy(src->xp, dst->xp);
    fp2copy(src->xq, dst->xq);
    fp2copy(src->xpq, dst->xpq);
}

static void swap_CurveAndPointsSIDH(CurveAndPointsSIDH *curveA, CurveAndPointsSIDH *curveB)
{
    CurveAndPointsSIDH *tmp = malloc(sizeof(CurveAndPointsSIDH));
    copy_CurveAndPointsSIDH(tmp, curveA);    // tmp = curveA
    copy_CurveAndPointsSIDH(curveA, curveB); // curveA = curveB
    copy_CurveAndPointsSIDH(curveB, tmp);    // curveB = tmp
    free(tmp);
}

static void precomputation_kernel(CurveAndPointsSIDH* list_curves, CurveAndPointsSIDH current_curve, size_t node_num, size_t h, unsigned long max_depth, unsigned long e)
{
    unsigned char k[4] = {4, 0, 0, 0}, l[4] = {0, 0, 0, 0};
    point_proj_t PpQ, A24, XP, XQ, XPQ, XR;
    point_proj_t A242, XP2, XQ2, XQ2_, XPQ2;

    // mark the input node as visited
    // calculate the node number on the tree based on inputs
    fp2copy(current_curve.a24, A24->X);
    fpcopy((digit_t *)&Montgomery_one, A24->Z[0]);
    fpzero(A24->Z[1]);
    fp2copy(current_curve.xp, XP->X);
    fpcopy((digit_t *)&Montgomery_one, XP->Z[0]); // P
    fpzero(XP->Z[1]);
    fp2copy(current_curve.xq, XQ->X);
    fpcopy((digit_t *)&Montgomery_one, XQ->Z[0]);
    fpzero(XQ->Z[1]);
    fp2copy(current_curve.xpq, XPQ->X);
    fpcopy((digit_t *)&Montgomery_one, XPQ->Z[0]); // P-Q
    fpzero(XPQ->Z[1]);

    fp2copy(XQ->X, XQ2->X);
    fpcopy((digit_t *)&Montgomery_one, XQ2->Z[0]);
    fpzero(XQ2->Z[1]);

    fp2copy(XP->X, PpQ->X);
    fpcopy((digit_t *)&Montgomery_one, PpQ->Z[0]);
    fpzero(PpQ->Z[1]);

    // Sum PpQ = x(P + Q)
    xDBLADD(XQ2, PpQ, XPQ->X, A24->X);
    fp2inv_mont(PpQ->Z);
    fp2mul_mont(PpQ->X, PpQ->Z, PpQ->X);

    xDBL_affine(XQ2, XQ2, A24->X);
    // generate more tasks that go deeper down the tree
    for (size_t i = 0; i < BRANCHING_FACTOR; i++)
    {
        size_t children_node_nums = node_num * BRANCHING_FACTOR + i;

        CurveAndPointsSIDH new_child;
        k[0] = (unsigned char)i;
        l[0] = (unsigned char)(4 - i);

        LadderThreePtSIKE(XP2, A24->X, XP->X, XQ->X, XPQ->X, k, 2);
        LadderThreePtSIKE(XPQ2, A24->X, XP->X, XQ->X, PpQ->X, l, 3);
        xDBLe_affine(XP2, XR, A24->X, e - 2 - 2 - 2 * h);

        IsogenyWithPoints(A242, XP2, XQ2_, XPQ2, A24->X, XR, XQ2, NULL, 0);

        if (h != (e - 4) / 2)
        {
            FourwayInv(A242->Z, XP2->Z, XQ2_->Z, XPQ2->Z);
            fp2mul_mont(A242->X, A242->Z, new_child.a24);
            fp2mul_mont(XP2->X, XP2->Z, new_child.xp);
            fp2mul_mont(XQ2_->X, XQ2_->Z, new_child.xq);
            fp2mul_mont(XPQ2->X, XPQ2->Z, new_child.xpq);
        }
        else
        {
            fp2add(A242->X, A242->X, A242->X);
            fp2sub(A242->X, A242->Z, A242->X);
            fp2add(A242->X, A242->X, A242->X);
            j_inv(A242->X, A242->Z, new_child.a24);

            /* Frobenius */
            fp2correction(new_child.a24);
            // if ((list_curves[index_w].a24)[1][0] & 1)
            //     fpneg((list_curves[index_w].a24)[1]);
        }

        if (h < max_depth - 1)
        {
#pragma omp task
            {
                precomputation_kernel(list_curves, new_child, children_node_nums, h + 1, max_depth, e);
            }
        } else {
			size_t index_w = children_node_nums;
			list_curves[index_w] = new_child;
        }
    }
}

void PrecompRightCurve_TaskBased_secondVer(CurveAndPointsSIDH *E1, CurveAndPointsSIDH *RightE, unsigned long delta, unsigned long e)
{
    /* First 1 step with kernel x = 1 */
	CurveAndPointsSIDH EA_R;

    point_proj_t PpQ, A24, XP, XQ, XPQ;
    unsigned char k[4] = {4, 0, 0, 0};

    fp2copy(RightE->xp, XP->X);
    fpcopy((digit_t *)&Montgomery_one, XP->Z[0]);
    fpzero(XP->Z[1]);

    fp2copy(RightE->xq, PpQ->X);
    fpcopy((digit_t *)&Montgomery_one, PpQ->Z[0]);
    fpzero(PpQ->Z[1]);

    xDBLADD(XP, PpQ, RightE->xpq, RightE->a24); // 2P, P+Q
    fp2inv_mont(PpQ->Z);
    fp2mul_mont(PpQ->X, PpQ->Z, PpQ->X); // P+Q
    xDBL_affine(XP, XP, RightE->a24);    // 4P

    LadderThreePtSIKE(XPQ, RightE->a24, RightE->xq, RightE->xp, PpQ->X, k, 3);

    fp2copy(RightE->xq, XQ->X);
    fpcopy((digit_t *)&Montgomery_one, XQ->Z[0]);
    fpzero(XQ->Z[1]);

    GetFourIsogenyWithKernelXeqZ(A24, RightE->a24);
    EvalFourIsogenyWithKernelXeqZ(XP, RightE->a24);
    EvalFourIsogenyWithKernelXeqZ(XQ, RightE->a24);
    EvalFourIsogenyWithKernelXeqZ(XPQ, RightE->a24);

    FourwayInv(A24->Z, XP->Z, XQ->Z, XPQ->Z);
    fp2mul_mont(A24->X, A24->Z, EA_R.a24);
    fp2mul_mont(XP->X, XP->Z, EA_R.xq);
    fp2mul_mont(XQ->X, XQ->Z, EA_R.xp); /* Swapping the values of P and Q */
    fp2mul_mont(XPQ->X, XPQ->Z, EA_R.xpq);

#pragma omp parallel
    {
#pragma omp single
        {
#pragma omp task
            {
                precomputation_kernel(E1, EA_R, 0, 0, delta / 2, e);
            }
        }
    }

#ifdef PRECOMP_OUTPUT_TO_FILE
#ifdef BIG_MEM_FOR_OUTPUT
    FILE *fpIn = fopen("PrecompRightCurve_TaskBased_secondVer.txt", "w");
#else
    FILE *fpIn = fopen("PrecompRightCurve_TaskBased_secondVer.txt", "w");
#endif
    for (size_t i = 0; i < (size_t)pow(2, delta); i++)
    {
        fprintf(fpIn, "%" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " \n", E1[i].a24[0][0], E1[i].a24[0][1], E1[i].a24[1][0], E1[i].a24[1][1]);
    }
#endif

#ifdef PERMUTATION
    // permutate
    for (size_t i = 0; i < (size_t)pow(2, delta); i++)
    {
        swap_CurveAndPointsSIDH(&E1[i], &E1[permutation_delta[i]]);
    }
#endif
//    free(wholeTree);
#ifdef PRECOMP_OUTPUT_TO_FILE
    fclose(fpIn);
#endif
}

void PrecompLeftCurve_TaskBased_secondVer(CurveAndPointsSIDH *E0, CurveAndPointsSIDH *LeftE, unsigned long delta, unsigned long e)
{
	CurveAndPointsSIDH E_left, E_right;
	CurveAndPointsSIDH *right_tree = &E0[(size_t)pow(2, delta)];
	CurveAndPointsSIDH* left_tree = E0;

    unsigned long h, j, index_r, index_w;
    unsigned char k[4], l[4]; /* Some max length.... should not be more */
    point_proj_t PpQ, A24, XP, XQ, XPQ, XR;
    point_proj_t A242, XP2, XQ2, XQ2_, XPQ2;

    /* First step bit = 0 */
    fp2copy(LeftE->xp, XP->X);
    fp2copy(LeftE->xq, XQ->X);
    fp2copy(LeftE->xpq, XPQ->X);

    fpcopy((digit_t *)&Montgomery_one, XP->Z[0]);
    fpzero(XP->Z[1]);
    fpcopy((digit_t *)&Montgomery_one, XQ->Z[0]);
    fpzero(XQ->Z[1]);
    fpcopy((digit_t *)&Montgomery_one, XPQ->Z[0]);
    fpzero(XPQ->Z[1]);

    xDBLADD(XQ, XPQ, LeftE->xp, LeftE->a24);
    xDBLe_affine(XP, XR, LeftE->a24, e - 2);

    GetTwoIsogenyWithXneZ(XR, A24);
    EvalTwoIsogenyWithXneZ(XR, XP);
    EvalTwoIsogenyWithXneZ(XR, XQ);
    EvalTwoIsogenyWithXneZ(XR, XPQ);

    FourwayInv(A24->Z, XP->Z, XQ->Z, XPQ->Z);
    fp2mul_mont(A24->X, A24->Z, E_left.a24);
    fp2mul_mont(XP->X, XP->Z, E_left.xp);
    fp2mul_mont(XQ->X, XQ->Z, E_left.xq);
    fp2mul_mont(XPQ->X, XPQ->Z, E_left.xpq);

    /* First step bit = 1 */
    /* Notice almost the same as above.. could share code */
    fp2copy(LeftE->xp, XP->X);
    fp2copy(LeftE->xq, XQ->X);
    fp2copy(LeftE->xpq, XPQ->X);

    fpcopy((digit_t *)&Montgomery_one, XP->Z[0]);
    fpzero(XP->Z[1]);
    fpcopy((digit_t *)&Montgomery_one, XQ->Z[0]);
    fpzero(XQ->Z[1]);
    fpcopy((digit_t *)&Montgomery_one, XPQ->Z[0]);
    fpzero(XPQ->Z[1]);

    xDBLADD(XQ, XP, LeftE->xpq, LeftE->a24);
    xDBLe_affine(XP, XR, LeftE->a24, e - 2);

    GetTwoIsogenyWithXneZ(XR, A24);
    EvalTwoIsogenyWithXneZ(XR, XP);
    EvalTwoIsogenyWithXneZ(XR, XQ);
    EvalTwoIsogenyWithXneZ(XR, XPQ);

    FourwayInv(A24->Z, XP->Z, XQ->Z, XPQ->Z);
    fp2mul_mont(A24->X, A24->Z, E_right.a24);
    fp2mul_mont(XP->X, XP->Z, E_right.xp);
    fp2mul_mont(XQ->X, XQ->Z, E_right.xq);
    fp2mul_mont(XPQ->X, XPQ->Z, E_right.xpq);

#pragma omp parallel
    {
#pragma omp single
        {
#pragma omp task
            {
                precomputation_kernel(left_tree, E_left, 0, 0, delta / 2, e);
            }
#pragma omp task
            {
                precomputation_kernel(right_tree, E_right, 0, 0, delta / 2, e);
            }
        }
    }

#ifdef PRECOMP_OUTPUT_TO_FILE
#ifdef BIG_MEM_FOR_OUTPUT
    FILE *fpIn = fopen("PrecompLeftCurve_TaskBased_secondVer.txt", "w");
#else
    FILE *fpIn = fopen("PrecompLeftCurve_TaskBased_secondVer.txt", "w");
#endif
    for (size_t i = 0; i < (size_t)pow(2, delta); i++)
    {
#pragma omp critical
        {
            fprintf(fpIn, "%" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " ", right_tree[i].a24[0][0], right_tree[i].a24[0][1], right_tree[i].a24[1][0], right_tree[i].a24[1][1]);
        }
    }

    for (size_t i = 0; i < (size_t)pow(2, delta); i++)
    {
#pragma omp critical
        {
            fprintf(fpIn, "%" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " ", E0[i].a24[0][0], E0[i].a24[0][1], E0[i].a24[1][0], E0[i].a24[1][1]);
        }
    }
#endif

#ifdef PERMUTATION
    // permutate
    for (size_t i = 0; i < (size_t)pow(2, delta); i++)
    {
        swap_CurveAndPointsSIDH(&E0[i], &E0[permutation_delta[i]]);
    }

    for (size_t i = 0; i < (size_t)pow(2, delta); i++)
    {
        swap_CurveAndPointsSIDH(&E0[i + (size_t)pow(2, delta)], &E0[permutation_delta[i] + (size_t)pow(2, delta)]);
    }
#endif
#ifdef PRECOMP_OUTPUT_TO_FILE
    fclose(fpIn);
#endif
}

void PrecompRightCurve(CurveAndPointsSIDH *E1, CurveAndPointsSIDH *RightE, unsigned long delta, unsigned long e)
{
    printf("Right curve\n");
#ifdef PRINT_CURVE_INFO
    print_curve_parameter(*RightE);
#endif
    /* First 1 step with kernel x = 1 */
    unsigned long h, j, index_r, index_w;
    unsigned char k[4] = {4, 0, 0, 0}, l[4] = {0, 0, 0, 0};
    point_proj_t PpQ, A24, XP, XQ, XPQ, XR;
    point_proj_t A242, XP2, XQ2, XQ2_, XPQ2;

    fp2copy(RightE->xp, XP->X);
    fpcopy((digit_t *)&Montgomery_one, XP->Z[0]);
    fpzero(XP->Z[1]);

    fp2copy(RightE->xq, PpQ->X);
    fpcopy((digit_t *)&Montgomery_one, PpQ->Z[0]);
    fpzero(PpQ->Z[1]);

    xDBLADD(XP, PpQ, RightE->xpq, RightE->a24); // 2P, P+Q
    fp2inv_mont(PpQ->Z);
    fp2mul_mont(PpQ->X, PpQ->Z, PpQ->X); // P+Q
    xDBL_affine(XP, XP, RightE->a24);    // 4P

    LadderThreePtSIKE(XPQ, RightE->a24, RightE->xq, RightE->xp, PpQ->X, k, 3);

    fp2copy(RightE->xq, XQ->X);
    fpcopy((digit_t *)&Montgomery_one, XQ->Z[0]);
    fpzero(XQ->Z[1]);

    GetFourIsogenyWithKernelXeqZ(A24, RightE->a24);
    EvalFourIsogenyWithKernelXeqZ(XP, RightE->a24);
    EvalFourIsogenyWithKernelXeqZ(XQ, RightE->a24);
    EvalFourIsogenyWithKernelXeqZ(XPQ, RightE->a24);

    FourwayInv(A24->Z, XP->Z, XQ->Z, XPQ->Z);
    fp2mul_mont(A24->X, A24->Z, E1[0].a24);
    fp2mul_mont(XP->X, XP->Z, E1[0].xq);
    fp2mul_mont(XQ->X, XQ->Z, E1[0].xp); /* Swapping the values of P and Q */
    fp2mul_mont(XPQ->X, XPQ->Z, E1[0].xpq);

#ifdef PRINT_CURVE_INFO
    printf("EA/R:\n");
    print_curve_parameter(E1[0]);
    return;
#endif
/* Now delta sized pre-computation */
#ifdef PRECOMP_OUTPUT_TO_FILE
#ifdef BIG_MEM_FOR_OUTPUT
    FILE *fpIn = fopen("PrecompRightCurve.txt", "w");
#else
    FILE *fpIn = fopen("PrecompRightCurve.txt", "w");
#endif
#endif
    for (h = 0; h < delta / 2; h++)
    {
        for (index_r = 0; index_r < (unsigned long)(1 << 2 * h); index_r++)
        {
            // Get curve and points P, Q, P-Q at index_r
            fp2copy(E1[index_r].a24, A24->X);
            fpcopy((digit_t *)&Montgomery_one, A24->Z[0]);
            fpzero(A24->Z[1]);
            fp2copy(E1[index_r].xp, XP->X);
            fpcopy((digit_t *)&Montgomery_one, XP->Z[0]); // P
            fpzero(XP->Z[1]);
            fp2copy(E1[index_r].xq, XQ->X);
            fpcopy((digit_t *)&Montgomery_one, XQ->Z[0]);
            fpzero(XQ->Z[1]);
            fp2copy(E1[index_r].xpq, XPQ->X);
            fpcopy((digit_t *)&Montgomery_one, XPQ->Z[0]); // P-Q
            fpzero(XPQ->Z[1]);

            fp2copy(XQ->X, XQ2->X);
            fpcopy((digit_t *)&Montgomery_one, XQ2->Z[0]);
            fpzero(XQ2->Z[1]);

            fp2copy(XP->X, PpQ->X);
            fpcopy((digit_t *)&Montgomery_one, PpQ->Z[0]);
            fpzero(PpQ->Z[1]);

            // Sum PpQ = x(P + Q)
            xDBLADD(XQ2, PpQ, XPQ->X, A24->X);
            fp2inv_mont(PpQ->Z);
            fp2mul_mont(PpQ->X, PpQ->Z, PpQ->X);

            xDBL_affine(XQ2, XQ2, A24->X);
#ifdef PRINT_CURVE_INFO
            // if (index_r == 2 && h == 5)
            // {
            // printf("index_r: %ld\n", index_r);
            // print_curve_parameter(E1[index_r]);
            // }
#endif
            for (j = 0; j < 4; j++)
            {
                k[0] = (unsigned char)j;
                l[0] = (unsigned char)(4 - j);

                LadderThreePtSIKE(XP2, A24->X, XP->X, XQ->X, XPQ->X, k, 2);
                LadderThreePtSIKE(XPQ2, A24->X, XP->X, XQ->X, PpQ->X, l, 3);
                xDBLe_affine(XP2, XR, A24->X, e - 2 - 2 - 2 * h);

                IsogenyWithPoints(A242, XP2, XQ2_, XPQ2, A24->X, XR, XQ2, NULL, 0);

                index_w = index_r + j * (1 << (2 * h));

                // #ifdef PRINT_CURVE_INFO
                // if (h == delta / 2 - 1)
                // {
                // Write curve and points at index_w
                if (2 * h + 2 != e - 2)
                {
                    FourwayInv(A242->Z, XP2->Z, XQ2_->Z, XPQ2->Z);
                    fp2mul_mont(A242->X, A242->Z, E1[index_w].a24);
                    fp2mul_mont(XP2->X, XP2->Z, E1[index_w].xp);
                    fp2mul_mont(XQ2_->X, XQ2_->Z, E1[index_w].xq);
                    fp2mul_mont(XPQ2->X, XPQ2->Z, E1[index_w].xpq);
#ifdef PRINT_CURVE_INFO
                    // if (h == delta / 2 - 1)
                    // {
                    // printf("\th = %ld index_r = %ld index_w = %ld\n", h, index_r, index_w);

                    // f2elm_t a24;
                    // printf("index_r: %ld\n", index_r);
                    // print_curve_parameter(E1[index_r]);

                    // printf("index_r: %ld index_w: %ld\n", index_r, index_w);
                    // print_curve_parameter(E1[index_w]);
                    // printf("kernel x-coord:\n");
                    // f2elm_t zR_inv;
                    // f2elm_t xR;
                    // fp2copy(XR[0].Z, zR_inv);
                    // fp2inv_mont(zR_inv);
                    // fp2mul_mont(XR[0].X, zR_inv, xR);
                    // from_mont(xR[0], a24[0]);
                    // from_mont(xR[1], a24[1]);
                    // print_f2elm(a24);
                    // if (index_w == 1026 && h == 5)
                    // {
                    //     getchar();
                    // }
                    // }
#endif
                }
                else
                {
                    //     /* Just store j-invariants to make it faster */
                    fp2add(A242->X, A242->X, A242->X);
                    fp2sub(A242->X, A242->Z, A242->X);
                    fp2add(A242->X, A242->X, A242->X);

                    // fp2inv_mont(A242->Z);
                    // fp2mul_mont(A242->X, A242->Z, E1[index_w].a24);
                    j_inv(A242->X, A242->Z, E1[index_w].a24);

                    /* Frobenius */
                    fp2correction(E1[index_w].a24);
                    // if ((E1[index_w].a24)[1][0] & 1)
                    //     fpneg((E1[index_w].a24)[1]);
#ifdef PRINT_CURVE_INFO
                    // printf("index_r: %ld leaf-level index_w: %ld\n", index_r, index_w);
                    // print_curve_parameter(E1[index_w]);
#endif
                }

                // if (index_w > max_idx_w)
                // {
                //     max_idx_w = index_w;
                //     printf("\thhhh\n");
                // }
                if (h == delta / 2 - 1)
                {
                    //                     f2elm_t a24;
                    //                     from_mont(E1[index_w].a24[0], a24[0]);
                    //                     from_mont(E1[index_w].a24[1], a24[1]);
#ifdef PRECOMP_OUTPUT_TO_FILE
                    fprintf(fpIn, "%" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " \n", E1[index_w].a24[0][0], E1[index_w].a24[0][1], E1[index_w].a24[1][0], E1[index_w].a24[1][1]);
                    // fprintf(fpIn, "%" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " \n", E1[index_w].xp[0][0], E1[index_w].xp[0][1], E1[index_w].xp[1][0], E1[index_w].xp[1][1]);
#endif

                    //                     // f2elm_t a24_to_search_for = {0x3e2d8532aa17fa6, 0x6, 0x47c7968f1bd62ce5, 0xb};
                    //                     // if (a24[0][0] == a24_to_search_for[0][0] &&
                    //                     //     a24[0][1] == a24_to_search_for[0][1] &&
                    //                     //     a24[1][0] == a24_to_search_for[1][0] &&
                    //                     //     a24[1][1] == a24_to_search_for[1][1])
                    //                     // {
                    //                     //     printf("\tfound\n");
                    //                     //     printf("\th = %ld\n", h);
                    //                     //     printf("\t\tindex_w = %ld\n", index_w);
                    //                     //     print_curve_parameter(E1[index_w]);
                    //                     //     printf("\t\tindex_r = %ld\n", index_r);
                    //                     //     print_curve_parameter(E1[index_r]);
                    //                     // }
                }
            }
        }
    }
#ifdef PRECOMP_OUTPUT_TO_FILE
    fclose(fpIn);
#endif

#ifdef PERMUTATION
    // permutate
    for (size_t i = 0; i < (size_t)pow(2, delta); i++)
    {
        swap_CurveAndPointsSIDH(&E1[i], &E1[permutation_delta[i]]);
    }
#endif
}

void PrecompLeftCurve(CurveAndPointsSIDH *E0, CurveAndPointsSIDH *LeftE, unsigned long delta, unsigned long e)
{
    printf("Left curve\n");
#ifdef PRINT_CURVE_INFO
    print_curve_parameter(*LeftE);
#endif
    unsigned long h, j, index_r, index_w;
    unsigned char k[4], l[4]; /* Some max length.... should not be more */
    point_proj_t PpQ, A24, XP, XQ, XPQ, XR;
    point_proj_t A242, XP2, XQ2, XQ2_, XPQ2;

    /* First step bit = 0 */
    fp2copy(LeftE->xp, XP->X);
    fp2copy(LeftE->xq, XQ->X);
    fp2copy(LeftE->xpq, XPQ->X);

    fpcopy((digit_t *)&Montgomery_one, XP->Z[0]);
    fpzero(XP->Z[1]);
    fpcopy((digit_t *)&Montgomery_one, XQ->Z[0]);
    fpzero(XQ->Z[1]);
    fpcopy((digit_t *)&Montgomery_one, XPQ->Z[0]);
    fpzero(XPQ->Z[1]);

    xDBLADD(XQ, XPQ, LeftE->xp, LeftE->a24);
    xDBLe_affine(XP, XR, LeftE->a24, e - 2);

    GetTwoIsogenyWithXneZ(XR, A24);
    EvalTwoIsogenyWithXneZ(XR, XP);
    EvalTwoIsogenyWithXneZ(XR, XQ);
    EvalTwoIsogenyWithXneZ(XR, XPQ);

    FourwayInv(A24->Z, XP->Z, XQ->Z, XPQ->Z);
    fp2mul_mont(A24->X, A24->Z, E0[0].a24);
    fp2mul_mont(XP->X, XP->Z, E0[0].xp);
    fp2mul_mont(XQ->X, XQ->Z, E0[0].xq);
    fp2mul_mont(XPQ->X, XPQ->Z, E0[0].xpq);
#ifdef PRINT_CURVE_INFO
    printf("index_r: E0[0]\n");
    print_curve_parameter(E0[0]);
#endif
    /* First step bit = 1 */
    /* Notice almost the same as above.. could share code */
    fp2copy(LeftE->xp, XP->X);
    fp2copy(LeftE->xq, XQ->X);
    fp2copy(LeftE->xpq, XPQ->X);

    fpcopy((digit_t *)&Montgomery_one, XP->Z[0]);
    fpzero(XP->Z[1]);
    fpcopy((digit_t *)&Montgomery_one, XQ->Z[0]);
    fpzero(XQ->Z[1]);
    fpcopy((digit_t *)&Montgomery_one, XPQ->Z[0]);
    fpzero(XPQ->Z[1]);

    xDBLADD(XQ, XP, LeftE->xpq, LeftE->a24);
    xDBLe_affine(XP, XR, LeftE->a24, e - 2);

    GetTwoIsogenyWithXneZ(XR, A24);
    EvalTwoIsogenyWithXneZ(XR, XP);
    EvalTwoIsogenyWithXneZ(XR, XQ);
    EvalTwoIsogenyWithXneZ(XR, XPQ);

    FourwayInv(A24->Z, XP->Z, XQ->Z, XPQ->Z);
    fp2mul_mont(A24->X, A24->Z, E0[1].a24);
    fp2mul_mont(XP->X, XP->Z, E0[1].xp);
    fp2mul_mont(XQ->X, XQ->Z, E0[1].xq);
    fp2mul_mont(XPQ->X, XPQ->Z, E0[1].xpq);
#ifdef PRINT_CURVE_INFO
    printf("index_r: E0[1]\n");
    print_curve_parameter(E0[1]);
    return;
#endif
#ifdef PRECOMP_OUTPUT_TO_FILE
#ifdef BIG_MEM_FOR_OUTPUT
    FILE *fpIn = fopen("PrecompLeftCurve.txt", "w");
#else
    FILE *fpIn = fopen("PrecompLeftCurve.txt", "w");
#endif
#endif
    if (delta != 0)
    {
        for (h = 0; h < delta / 2; h++)
        {
            for (index_r = 0; index_r < (unsigned long)(2 * (1 << 2 * h)); index_r++)
            {
                // Get curve and points P, Q, P-Q at index_r
                fp2copy(E0[index_r].a24, A24->X);
                fpcopy((digit_t *)&Montgomery_one, A24->Z[0]);
                fpzero(A24->Z[1]);
                fp2copy(E0[index_r].xp, XP->X);
                fpcopy((digit_t *)&Montgomery_one, XP->Z[0]); // P
                fpzero(XP->Z[1]);
                fp2copy(E0[index_r].xq, XQ->X);
                fpcopy((digit_t *)&Montgomery_one, XQ->Z[0]);
                fpzero(XQ->Z[1]);
                fp2copy(E0[index_r].xpq, XPQ->X);
                fpcopy((digit_t *)&Montgomery_one, XPQ->Z[0]); // P-Q
                fpzero(XPQ->Z[1]);

                fp2copy(XQ->X, XQ2->X);
                fpcopy((digit_t *)&Montgomery_one, XQ2->Z[0]);
                fpzero(XQ2->Z[1]);

                fp2copy(XP->X, PpQ->X);
                fpcopy((digit_t *)&Montgomery_one, PpQ->Z[0]);
                fpzero(PpQ->Z[1]);

                // Sum PpQ = x(P + Q)
                xDBLADD(XQ2, PpQ, XPQ->X, A24->X);
                fp2inv_mont(PpQ->Z);
                fp2mul_mont(PpQ->X, PpQ->Z, PpQ->X);

                xDBL_affine(XQ2, XQ2, A24->X);

#ifdef PRINT_CURVE_INFO
                // bool print_index_r = false;
                // if (!print_index_r)
                // {
                //     printf("index_r: %ld\n", index_r);
                //     print_curve_parameter(E0[index_r]);
                //     print_index_r = true;
                // }
#endif
                for (j = 0; j < 4; j++)
                {
                    k[0] = (unsigned char)j;
                    l[0] = (unsigned char)(4 - j);

                    LadderThreePtSIKE(XP2, A24->X, XP->X, XQ->X, XPQ->X, k, 2);
                    LadderThreePtSIKE(XPQ2, A24->X, XP->X, XQ->X, PpQ->X, l, 3);
                    xDBLe_affine(XP2, XR, A24->X, e - 2 - 2 - 2 * h);

                    IsogenyWithPoints(A242, XP2, XQ2_, XPQ2, A24->X, XR, XQ2, NULL, 0);
#ifdef PRINT_CURVE_INFO
                    // if (h == delta / 2 - 1)
                    // {
                    //     //     printf("index_r: %ld\n", index_r);
                    //     //     print_curve_parameter(E0[index_r]);

                    //     f2elm_t a24;
                    //     printf("kernel x-coord:\n");
                    //     f2elm_t zR_inv;
                    //     f2elm_t xR;
                    //     fp2copy(XR[0].Z, zR_inv);
                    //     fp2inv_mont(zR_inv);
                    //     fp2mul_mont(XR[0].X, zR_inv, xR);
                    //     from_mont(xR[0], a24[0]);
                    //     from_mont(xR[1], a24[1]);
                    //     print_f2elm(a24);
                    // }
#endif
                    index_w = index_r + j * 2 * (1 << (2 * h));
                    // write curve and points at index_w
                    if (2 * h + 2 != e - 2)
                    {
                        FourwayInv(A242->Z, XP2->Z, XQ2_->Z, XPQ2->Z);
                        fp2mul_mont(A242->X, A242->Z, E0[index_w].a24);
                        fp2mul_mont(XP2->X, XP2->Z, E0[index_w].xp);
                        fp2mul_mont(XQ2_->X, XQ2_->Z, E0[index_w].xq);
                        fp2mul_mont(XPQ2->X, XPQ2->Z, E0[index_w].xpq);

#ifdef PRINT_CURVE_INFO
                        // if (h == delta / 2 - 1)
                        // {
                        // printf("\th = %ld index_r = %ld index_w = %ld\n", h, index_r, index_w);

                        // f2elm_t a24;
                        // printf("index_r: %ld\n", index_r);
                        // print_curve_parameter(E0[index_r]);

                        // printf("index_r: %ld index_w: %ld\n", index_r, index_w);
                        // printf("index_w: %ld\n", index_w);
                        // print_curve_parameter(E0[index_w]);
                        // printf("kernel x-coord:\n");
                        // f2elm_t zR_inv;
                        // f2elm_t xR;
                        // fp2copy(XR[0].Z, zR_inv);
                        // fp2inv_mont(zR_inv);
                        // fp2mul_mont(XR[0].X, zR_inv, xR);
                        // from_mont(xR[0], a24[0]);
                        // from_mont(xR[1], a24[1]);
                        // print_f2elm(a24);
                        // if (index_w == 1026 && h == 5)
                        // {
                        //     getchar();
                        // }
                        // }
#endif
                    }
                    else
                    {
                        /* Just store j-invariants to make it faster */
                        fp2add(A242->X, A242->X, A242->X);
                        fp2sub(A242->X, A242->Z, A242->X);
                        fp2add(A242->X, A242->X, A242->X);
                        j_inv(A242->X, A242->Z, E0[index_w].a24);

                        // fp2inv_mont(A242->Z);
                        // fp2mul_mont(A242->X, A242->Z, E0[index_w].a24);
#ifdef PRINT_CURVE_INFO
                        // printf("index_r: %ld leaf-level index_w: %ld before correction\n", index_r, index_w);
                        // print_curve_parameter(E0[index_w]);
#endif
                        /* Frobenius */
                        fp2correction(E0[index_w].a24);
#ifdef PRINT_CURVE_INFO
                        // printf("index_r: %ld leaf-level index_w: %ld after correction\n", index_r, index_w);
                        // print_curve_parameter(E0[index_w]);
                        // if (index_w == 8193)
                        // {
                        //     printf("test0:\n");
                        //     print_felm_hex((E0[index_w].a24)[0], false);
                        //     print_felm_hex((E0[index_w].a24)[1], true);
                        //     getchar();
                        // }
#endif
                        // if ((E0[index_w].a24)[1][0] & 1)
                        // {
                        //     // printf("test1:\n");
                        //     // print_felm_hex((E0[index_w].a24)[0], false);
                        //     // print_felm_hex((E0[index_w].a24)[1], true);
                        //     // getchar();
                        //     fpneg((E0[index_w].a24)[1]);
                        // }
#ifdef PRINT_CURVE_INFO
                        // printf("leaf-level index_w: %ld\n", index_w);
                        // print_curve_parameter(E0[index_w]);
                        // printf("index_r: %ld leaf-level index_w: %ld\n", index_r, index_w);
                        // print_curve_parameter(E0[index_w]);
#endif
                    }
                    if (h == delta / 2 - 1)
                    {
#ifdef PRECOMP_OUTPUT_TO_FILE
                        //                         printf("index_r: %ld leaf-level index_w: %ld after Frobenius\n", index_r, index_w);
                        //                         print_curve_parameter(E0[index_w]);
                        //                         // f2elm_t a24;
                        //                         // from_mont(E0[index_w].a24[0], a24[0]);
                        //                         // from_mont(E0[index_w].a24[1], a24[1]);
                        fprintf(fpIn, "%" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " ", E0[index_w].a24[0][0], E0[index_w].a24[0][1], E0[index_w].a24[1][0], E0[index_w].a24[1][1]);
#endif
                    }
                }
            }
        }
    }
#ifdef PERMUTATION
    // permutate
    for (size_t i = 0; i < (size_t)pow(2, delta); i++)
    {
        swap_CurveAndPointsSIDH(&E0[i], &E0[permutation_delta[i]]);
    }

    for (size_t i = 0; i < (size_t)pow(2, delta); i++)
    {
        swap_CurveAndPointsSIDH(&E0[i + (size_t)pow(2, delta)], &E0[permutation_delta[i] + (size_t)pow(2, delta)]);
    }
#endif
#ifdef PRECOMP_OUTPUT_TO_FILE
    fclose(fpIn);
#endif
}

static void input_points_to_mont()
{
    printf("input_points_to_mont\n");
    // input points in F_p2 displayed in hexadecimal
    // byte sequences are stored in little-endian

    // My instance p36_22
    f2elm_t a24 = {0x0000000000000002, 0x0000000000000000, 0x0, 0x0};
    f2elm_t xp = {0x9E378723AD64D83D, 0x0000000000000908, 0x4E4CAA3D14E98878, 0x0000000000000028};
    f2elm_t xq = {0xAE3E679AB8A2B807, 0x00000000000004E7, 0x52074933F24453B9, 0x0000000000000C3C};
    f2elm_t xpq = {0xD65F8A020C61E222, 0x00000000000006C4, 0x1C5C04FC4AF1FE97, 0x0000000000000ACA};

    f2elm_t a24_mont;
    f2elm_t xp_mont;
    f2elm_t xq_mont;
    f2elm_t xpq_mont;

    to_mont(a24[0], a24_mont[0]);
    to_mont(a24[1], a24_mont[1]);
    to_mont(xp[0], xp_mont[0]);
    to_mont(xp[1], xp_mont[1]);
    to_mont(xq[0], xq_mont[0]);
    to_mont(xq[1], xq_mont[1]);
    to_mont(xpq[0], xpq_mont[0]);
    to_mont(xpq[1], xpq_mont[1]);
    printf("E6:\n");
    printf(".a24=");
    print_f2elm(a24_mont);
    printf(".xp=");
    print_f2elm(xp_mont);
    printf(".xq=");
    print_f2elm(xq_mont);
    printf(".xpq=");
    print_f2elm(xpq_mont);

    f2elm_t EA_a24 = {0xFB300CA88F66649D, 0x0000000000000244, 0xC204A7131D24BE9D, 0x00000000000001EA};
    f2elm_t EA_xp = {0xA28C32A47C9ACB87, 0x00000000000001F5, 0xDE0E095912E4440F, 0x000000000000080A};
    f2elm_t EA_xq = {0x413F5FE8E2DB8A77, 0x000000000000097E, 0x3ECC91F46F1C8C90, 0x00000000000007F3};
    f2elm_t EA_xpq = {0x6E51BEA9542D99C7, 0x00000000000007C2, 0xFEB97159986468C8, 0x00000000000004BA};
    to_mont(EA_a24[0], a24_mont[0]);
    to_mont(EA_a24[1], a24_mont[1]);
    to_mont(EA_xp[0], xp_mont[0]);
    to_mont(EA_xp[1], xp_mont[1]);
    to_mont(EA_xq[0], xq_mont[0]);
    to_mont(EA_xq[1], xq_mont[1]);
    to_mont(EA_xpq[0], xpq_mont[0]);
    to_mont(EA_xpq[1], xpq_mont[1]);
    printf("EA:\n");
    printf(".a24=");
    print_f2elm(a24_mont);
    printf(".xp=");
    print_f2elm(xp_mont);
    printf(".xq=");
    print_f2elm(xq_mont);
    printf(".xpq=");
    print_f2elm(xpq_mont);
}

static void read_end_curve_parameter(instance_t *inst)
{
    f2elm_t a24;
    f2elm_t xp_standard;
    f2elm_t xq_standard;
    f2elm_t xpq_standard;
    f2elm_t j_inv_golden_coll;

    from_mont(inst->jinv[0], j_inv_golden_coll[0]);
    from_mont(inst->jinv[1], j_inv_golden_coll[1]);
    printf("j_inv:\n");
    print_f2elm(j_inv_golden_coll);

    printf("Left curve:\n");
    print_curve_parameter(inst->E[0]);

    printf("Right curve:\n");
    print_curve_parameter(inst->E[1]);
}

void init_shared_state_second_step(instance_t *inst, shared_state_t *S)
{
    unsigned long strat[250];
    unsigned int i;

    S->lenstrat = OptStrat(strat, (unsigned long)(inst->e - 2 - S->delta) / 2, 1, 1);
    S->strat = calloc(S->lenstrat, sizeof(unsigned long));
    for (i = 0; i < S->lenstrat; i++)
    {
        S->strat[i] = strat[i];
    }

    /* Initialize Hansel & Gretel */
    S->HANSEL_GRETEL = inst->HANSEL_GRETEL;
    S->MAX_CRUMBS = inst->MAX_CRUMBS;

    /* Initialize vOW params */
    fp2copy(inst->jinv, S->jinv); /* The solution */
    double THETA = inst->ALPHA * sqrt((double)S->MEMORY_SIZE / (double)((unsigned long)1 << (inst->e - 1)));
    S->MAX_STEPS = ceil(inst->GAMMA / THETA);
    S->MAX_DIST = (unsigned long)(inst->BETA * S->MEMORY_SIZE);
    S->MAX_FUNCTION_VERSIONS = 10000;
    S->DIST_BOUND = THETA * (1 << (S->NBITS_STATE - S->MEMORY_LOG_SIZE));
    // assumes we are not filling digit_t, should shift by min(sizeof(digit_t)*8, inst->NBITS_STATE) - S->MEMORY_LOG_SIZE
    assert(S->NBITS_STATE - S->MEMORY_LOG_SIZE > 0);

    /* Statistics */
    S->collect_vow_stats = false; // By default don't collect stats (=> terminate run when successful)
    initTree(&S->dist_cols);      // Initing even if not using
    S->success = false;
    S->wall_time = 0.;
    S->collisions = 0;
    S->mem_collisions = 0;
    S->dist_points = 0;
    S->number_steps_collect = 0;
    S->number_steps_locate = 0;
    S->number_steps = 0;
    S->initial_function_version = 1;
    S->final_avg_random_functions = 0.;

    /* resync */
    S->resync_frequency = 10;
    S->resync_cores = (uint8_t *)calloc(S->N_OF_CORES, sizeof(uint8_t));
}

void init_shared_state(instance_t *inst, shared_state_t *S)
{
    printf("e = %ld\n", inst->e);
    printf("Instance: %s\n", inst->MODULUS);

    /* Initialize state */
    S->instance = inst;
    S->NBITS_STATE = inst->e - 1;                                 /* Walk of size e from 2 sides*/
    S->NBYTES_STATE = ((S->NBITS_STATE + 7) / 8);                 /* Number of bytes needed for state */
    S->NWORDS_STATE = ((S->NBITS_STATE + RADIX64 - 1) / RADIX64); /* Number of words need for state, assumes 64-bit arch */
    S->NBITS_OVERFLOW = (S->NBITS_STATE % 8) == 0 ? 8 : (S->NBITS_STATE % 8);
    S->PRNG_SEED = (unsigned long)inst->PRNG_SEED;

    /* Initialize memory */
    S->MEMORY_LOG_SIZE = inst->MEMORY_LOG_SIZE;
    S->MEMORY_SIZE = (uint64_t)(1 << S->MEMORY_LOG_SIZE);
    S->MEMORY_SIZE_MASK = S->MEMORY_SIZE - 1;
    if (initialize_shared_memory(S) == EXIT_FAILURE)
    {
        printf("Error initialising shared memory\n");
        assert(0);
    }

    /* Initialize omp params */
    S->N_OF_CORES = inst->N_OF_CORES;

    /* Initialize pre-computed tables */
    S->delta = inst->delta;

    S->external_E[0] = false;
    S->external_E[1] = false;
    double start;
    S->E[0] = calloc(((size_t)1 << (S->delta + 1)) * 4 * 2 * (size_t)NWORDS_FIELD, sizeof(digit_t));
    S->E[1] = calloc(((size_t)1 << S->delta) * 4 * 2 * (size_t)NWORDS_FIELD, sizeof(digit_t));
    if (!S->E[0] || !S->E[1]){
        printf("calloc failed\n");
    }
}

void free_shared_state(shared_state_t *S)
{
    cleanup_shared_memory(S);
    free(S->strat);

    // Free internally allocated precomputation tables
    for (unsigned int i = 0; i < 2; i++)
    {
        if (!S->external_E[i])
        {
            free(S->E[i]);
        }
    }

    if (S->dist_cols.size != 0)
    {
        freeTree(S->dist_cols.root);
    }

    free(S->resync_cores);
}

void SampleSIDH(private_state_t *private_state)
{
    sample_prng(&private_state->prng_state, private_state->current.current_state.bytes, (unsigned long)private_state->NBYTES_STATE);

    private_state->current.current_steps = 0;
    fix_overflow(&private_state->current.current_state, private_state->NBYTES_STATE, private_state->NBITS_OVERFLOW);
    copy_st(&private_state->current.initial_state, &private_state->current.current_state, private_state->NWORDS_STATE);

    // Hansel & Gretel
    clean_private_state(private_state);
}

void LadderThreePtSIDH(point_proj_t R,
                       const f2elm_t a24,
                       const f2elm_t xp,
                       const f2elm_t xq,
                       const f2elm_t xpq,
                       const unsigned char c,
                       const unsigned char *m,
                       unsigned long nbits_state,
                       unsigned long delta)
{ // Non-constant time version that depends on size of k
    point_proj_t R0 = {0}, R2 = {0};
    unsigned int i, bit, msb;

    fp2copy(xq, R0->X);
    fp2copy(xp, R->X);
    fp2copy(xpq, R2->X);
    fpcopy((digit_t *)&Montgomery_one, (digit_t *)R0->Z);
    fpzero((digit_t *)(R0->Z)[1]);
    fpcopy((digit_t *)&Montgomery_one, (digit_t *)R->Z);
    fpzero((digit_t *)(R->Z)[1]);
    fpcopy((digit_t *)&Montgomery_one, (digit_t *)R2->Z);
    fpzero((digit_t *)(R2->Z)[1]);

    msb = GetMSBSIDH(m, nbits_state); /* Can skip top zeroes of k */
    for (i = 2 - c + delta; i < msb; i++)
    { // Ignore c
        bit = (m[i >> 3] >> (i & 0x07)) & 1;

        if (bit)
        {
            xDBLADD_SIDH(R0, R, R2->X, a24);
            fp2mul_mont(R->X, R2->Z, R->X);
        }
        else
        {
            xDBLADD_SIDH(R0, R2, R->X, a24);
            fp2mul_mont(R2->X, R->Z, R2->X);
        }
    }
}

void GetIsogeny(f2elm_t jinv,
                const f2elm_t a24,
                const f2elm_t xp,
                const f2elm_t xq,
                const f2elm_t xpq,
                const unsigned char c,
                const unsigned char *k,
                const unsigned long *strat,
                const unsigned long lenstrat,
                unsigned long nbits_state,
                unsigned long delta)
{
    point_proj_t R, A24, pts[MAX_INT_POINTS_ALICE];
    f2elm_t coeff[3];
    unsigned long i, row, index = 0, ii = 0, m, npts = 0, pts_index[MAX_INT_POINTS_ALICE];

    if (delta + 2 != nbits_state + 1)
    {
        /* Retrieve kernel point */
        LadderThreePtSIDH(R, a24, xp, xq, xpq, c, k, nbits_state, delta);

        fp2copy(a24, A24->X);
        fpcopy((digit_t *)&Montgomery_one, (digit_t *)A24->Z);
        fpzero((digit_t *)(A24->Z)[1]);

        /* NOTE: Make the first step use the fact that a24 is affine? */
        for (row = 1; row < lenstrat + 1; row++)
        {
            while (index < lenstrat + 1 - row)
            {
                fp2copy(R->X, pts[npts]->X);
                fp2copy(R->Z, pts[npts]->Z);
                pts_index[npts++] = index;
                m = strat[ii++];
                xDBLe(R, R, A24->X, A24->Z, (int)(2 * m));
                index += m;
            }
            GetFourIsogenyWithKernelXneZ(R, A24->X, A24->Z, coeff);
            for (i = 0; i < npts; i++)
            {
                EvalFourIsogenyWithKernelXneZ(pts[i], coeff);
            }

            fp2copy(pts[npts - 1]->X, R->X);
            fp2copy(pts[npts - 1]->Z, R->Z);
            index = pts_index[npts - 1];
            npts -= 1;
        }
        GetFourIsogenyWithKernelXneZ(R, A24->X, A24->Z, coeff);

        fp2add(A24->X, A24->X, A24->X);
        fp2sub(A24->X, A24->Z, A24->X);
        fp2add(A24->X, A24->X, A24->X);
        j_inv(A24->X, A24->Z, jinv);

        /* Frobenius */
        fp2correction(jinv);
        if (jinv[1][0] & 1)
        {
            fpneg(jinv[1]);
        }
    }
    else
    {
        fp2copy(a24, jinv);
    }
}

void UpdateStSIDH(unsigned char jinvariant[FP2_ENCODED_BYTES], st_t *r, const st_t *s, private_state_t *private_state)
{
    f2elm_t jinv;
    unsigned char c = GetC_SIDH(s);
    unsigned long index;

    /* Get the j-invariant of the corresponding curve */
    if (c == 0)
    {
        index = (s->words[0] >> 1) & ((1 << (private_state->delta + 1)) - 1);
        GetIsogeny(jinv, private_state->E[0][index].a24,
                   private_state->E[0][index].xp,
                   private_state->E[0][index].xq,
                   private_state->E[0][index].xpq,
                   c, s->bytes, private_state->strat, private_state->lenstrat, (unsigned long)private_state->NBITS_STATE, private_state->delta);
    }
    else
    {
        index = (s->words[0] >> 1) & ((1 << private_state->delta) - 1);
        GetIsogeny(jinv, private_state->E[1][index].a24,
                   private_state->E[1][index].xp,
                   private_state->E[1][index].xq,
                   private_state->E[1][index].xpq,
                   c, s->bytes, private_state->strat, private_state->lenstrat, (unsigned long)private_state->NBITS_STATE, private_state->delta);
    }

    /* Hash j into (c,b,k) */
    fp2_encode(jinv, jinvariant); /* Unique encoding (includes fpcorrection) */
    XOF(r->bytes, jinvariant, (unsigned long)private_state->NBYTES_STATE, FP2_ENCODED_BYTES, (unsigned long)private_state->function_version);
    fix_overflow(r, private_state->NBYTES_STATE, private_state->NBITS_OVERFLOW);
}
