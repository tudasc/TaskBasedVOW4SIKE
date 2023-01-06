#pragma once
void LadderThreePtSIDH(point_proj_t R,
                       const f2elm_t a24,
                       const f2elm_t xp,
                       const f2elm_t xq,
                       const f2elm_t xpq,
                       const unsigned char c,
                       const unsigned char *m,
                       unsigned long nbits_state,
                       unsigned long delta);
void LadderThreePtSIKE(point_proj_t R, const f2elm_t a24, const f2elm_t xp, const f2elm_t xq, const f2elm_t xpq, const unsigned char *m, unsigned long nbits_state);
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
                unsigned long delta);
void GetTwoIsogenyWithXneZ(const point_proj_t R, point_proj_t A24);
void EvalTwoIsogenyWithXneZ(const point_proj_t R, point_proj_t P);
void IsogenyWithPoints(point_proj_t A24, point_proj_t XP, point_proj_t XQ, point_proj_t XPQ,
                       const f2elm_t a24, const point_proj_t kernel_point, const point_proj_t XQinp,
                       const unsigned long *strat, const unsigned long lenstrat);

typedef struct shared_container
{
    CurveAndPointsSIDH *list_curves;
    uint64_t current_index;
    omp_lock_t my_lock;
} shared_container;

extern uint64_t *starting_index_for_each_level;
extern uint64_t max_index_w;
extern shared_container right_precomputation_tree;
extern shared_container E_left;
extern shared_container E_right;

void PrecompRightCurve(CurveAndPointsSIDH *E1, CurveAndPointsSIDH *RightE, unsigned long delta, unsigned long e);
void PrecompLeftCurve(CurveAndPointsSIDH *E0, CurveAndPointsSIDH *LeftE, unsigned long delta, unsigned long e);

void PrecompRightCurve_TaskBased_secondVer(CurveAndPointsSIDH *E1, CurveAndPointsSIDH *RightE, unsigned long delta, unsigned long e);
void PrecompLeftCurve_TaskBased_secondVer(CurveAndPointsSIDH *E0, CurveAndPointsSIDH *LeftE, unsigned long delta, unsigned long e);

void PrecompRightCurve_TaskBased_thirdVersion(CurveAndPointsSIDH *E1, CurveAndPointsSIDH *RightE, unsigned long delta, unsigned long e);
void PrecompLeftCurve_TaskBased_thirdVersion(CurveAndPointsSIDH *E0, CurveAndPointsSIDH *LeftE, unsigned long delta, unsigned long e);