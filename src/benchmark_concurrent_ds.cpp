#include "benchmark_concurrent_ds.h"
#include <omp.h>
#include <iostream>
#include <cassert>

using namespace std;

extern "C"
{
#include "P128_internal.h"
#include "sike_vow.h"
#include "sidh_vow_base.h"
#include "sike_vow_instances.h"
#include "memory.h"
#include "bintree.h"
}

#include <cstring>
#define BRANCHING_FACTOR 4

void precomputation_kernel_forth_version(CurveAndPointsSIDH current_curve, tbb_vector *container, size_t h, unsigned long max_depth, unsigned long e)
{
    unsigned char k[4] = {4, 0, 0, 0}, l[4] = {0, 0, 0, 0};
    point_proj_t PpQ, A24, XP, XQ, XPQ, XR;
    point_proj_t A242, XP2, XQ2, XQ2_, XPQ2;

    // mark the input node as visited
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
        CurveAndPointsSIDH new_child;
#ifdef DEBUG
#pragma omp critical
        printf("\t\tthread %d h = %ld read idx %ld store idx %ld\n", omp_get_thread_num(), h, index_r, index_w);
#endif
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
            // if ((new_child.a24)[1][0] & 1)
            //     fpneg((new_child.a24)[1]);
        }

        if (h < max_depth - 1)
        {
#pragma omp task
            {
                precomputation_kernel_forth_version(new_child, container, h + 1, max_depth, e);
            }
        }
        else
        {
            container->push_back(std::move(new_child));
        }
    }
}

void PrecompRightCurve_TaskBased_TBB_Version(CurveAndPointsSIDH *E1, CurveAndPointsSIDH *RightE, unsigned long delta, unsigned long e)
{
    tbb_vector temp_vect;
    temp_vect.reserve((size_t)pow(2,delta));

    CurveAndPointsSIDH current_curve;

    /* First 1 step with kernel x = 1 */
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
    fp2mul_mont(A24->X, A24->Z, current_curve.a24);
    fp2mul_mont(XP->X, XP->Z, current_curve.xq);
    fp2mul_mont(XQ->X, XQ->Z, current_curve.xp); /* Swapping the values of P and Q */
    fp2mul_mont(XPQ->X, XPQ->Z, current_curve.xpq);

#pragma omp parallel
    {
#pragma omp single
        {
#pragma omp task
            {
                precomputation_kernel_forth_version(current_curve, &temp_vect, 0, delta / 2, e);
            }
        }
    }

    size_t i = 0;
    for(const auto& it: temp_vect){
        memcpy(&E1[i], &it, 4 * 2 * (size_t)NWORDS_FIELD * sizeof(digit_t));
        ++i;
    }
    assert(i == (size_t)pow(2,delta));

#ifdef PRECOMP_OUTPUT_TO_FILE
#ifdef BIG_MEM_FOR_OUTPUT
    FILE *fpIn = fopen("PrecompRightCurve_TaskBased_TBB_Ver.txt", "w");
#else
    FILE *fpIn = fopen("PrecompRightCurve_TaskBased_TBB_Ver.txt", "w");
#endif
    for (size_t i = 0; i < (size_t)pow(2, delta); i++)
    {
        fprintf(fpIn, "%" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " ", E1[i].a24[0][0], E1[i].a24[0][1], E1[i].a24[1][0], E1[i].a24[1][1]);
    }
#endif
#ifdef PRECOMP_OUTPUT_TO_FILE
    fclose(fpIn);
#endif
}

void PrecompLeftCurve_TaskBased_TBB_Version(CurveAndPointsSIDH *E0, CurveAndPointsSIDH *LeftE, unsigned long delta, unsigned long e)
{
    tbb_vector left_sub_tree;
    left_sub_tree.reserve((size_t)pow(2,delta));
    tbb_vector right_sub_tree;
    right_sub_tree.reserve((size_t)pow(2,delta));

    CurveAndPointsSIDH current_curve_E_left;
    CurveAndPointsSIDH current_curve_E_right;

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
    fp2mul_mont(A24->X, A24->Z, current_curve_E_left.a24);
    fp2mul_mont(XP->X, XP->Z, current_curve_E_left.xp);
    fp2mul_mont(XQ->X, XQ->Z, current_curve_E_left.xq);
    fp2mul_mont(XPQ->X, XPQ->Z, current_curve_E_left.xpq);

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
    fp2mul_mont(A24->X, A24->Z, current_curve_E_right.a24);
    fp2mul_mont(XP->X, XP->Z, current_curve_E_right.xp);
    fp2mul_mont(XQ->X, XQ->Z, current_curve_E_right.xq);
    fp2mul_mont(XPQ->X, XPQ->Z, current_curve_E_right.xpq);
#pragma omp parallel
    {
#pragma omp single
        {
#pragma omp task
            {
                precomputation_kernel_forth_version(current_curve_E_left, &left_sub_tree, 0, delta / 2, e);
            }
#pragma omp task
            {
                precomputation_kernel_forth_version(current_curve_E_right, &right_sub_tree, 0, delta / 2, e);
            }
        }
    }
    // copy left part
    size_t i = 0;
    for (const auto& it : left_sub_tree)
    {
        memcpy(&E0[i], &it, 4 * 2 * (size_t)NWORDS_FIELD * sizeof(digit_t));
        ++i;
    }

    // copy right part
    for (const auto& it : right_sub_tree)
    {
        memcpy(&E0[i], &it, 4 * 2 * (size_t)NWORDS_FIELD * sizeof(digit_t));
        ++i;
    }
    assert(i == (size_t)pow(2,delta+1));
    CurveAndPointsSIDH* right_part_S_E0 = &E0[(size_t)pow(2,delta)];

#ifdef PRECOMP_OUTPUT_TO_FILE
#ifdef BIG_MEM_FOR_OUTPUT
    FILE *fpIn = fopen("PrecompLeftCurve_TaskBased_TBB_Ver.txt", "w");
#else
    FILE *fpIn = fopen("PrecompLeftCurve_TaskBased_TBB_Ver.txt", "w");
#endif
    for (size_t i = 0; i < (size_t)pow(2, delta); i++)
    {
        fprintf(fpIn, "%" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " ", right_part_S_E0[i].a24[0][0], right_part_S_E0[i].a24[0][1], right_part_S_E0[i].a24[1][0], right_part_S_E0[i].a24[1][1]);
    }

    for (size_t i = 0; i < (size_t)pow(2, delta); i++)
    {
        fprintf(fpIn, "%" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " ", E0[i].a24[0][0], E0[i].a24[0][1], E0[i].a24[1][0], E0[i].a24[1][1]);
    }
#endif
#ifdef PRECOMP_OUTPUT_TO_FILE
    fclose(fpIn);
#endif
}

void benchmark_precomputation(instance_t* inst, shared_state_t* S)
{
    double start;
    cout << "Precomputing...\n";
#ifdef TASK_BASED_PRECOMP_SECOND_VER
    start = omp_get_wtime();
    PrecompLeftCurve_TaskBased_secondVer(S->E[0], &inst->E[0], S->delta, inst->e);
    PrecompRightCurve_TaskBased_secondVer(S->E[1], &inst->E[1], S->delta, inst->e);
    printf("Task-based second-version precomputation: %f s\n", omp_get_wtime() - start);
#elif TASK_BASED_PRECOMP_TBB_VER
    start = omp_get_wtime();
    PrecompLeftCurve_TaskBased_TBB_Version(S->E[0], &inst->E[0], S->delta, inst->e);
    PrecompRightCurve_TaskBased_TBB_Version(S->E[1], &inst->E[1], S->delta, inst->e);
    printf("Task-based tbb-version precomputation: %f s\n", omp_get_wtime() - start);
#elif SERIAL_VER
    start = omp_get_wtime();
    PrecompRightCurve(S->E[1], &inst->E[1], S->delta, inst->e);
    PrecompLeftCurve(S->E[0], &inst->E[0], S->delta, inst->e);
    printf("Original precomputation: %f s\n", omp_get_wtime() - start);
#endif
    
    free(starting_index_for_each_level);

#ifdef MITM
    free(S->E[0]);
    free(S->E[1]);
    free(S);
    printf("Precomputation done!\n");
    return;
#endif
}
