#include <string.h>
#include <stdio.h>
#include "types/state.h"

#if defined(VOW_SIDH)
#elif defined(VOW_SIKE)
#else
#endif

#if defined(VOW_SIDH) || defined(VOW_SIKE)
void print_felm(felm_t *el)
{
    for (uint16_t i = 0; i < NWORDS_FIELD; i++)
        printf("%" PRIu64 " ", ((digit_t *)el)[i]);
    printf("\n");
}

void print_felm_hex(felm_t *el, bool img_part)
{
    if (!img_part)
    {
        for (uint16_t i = 0; i < NWORDS_FIELD; i++)
            printf("%#018lx, ", ((digit_t *)el)[i]);
    }
    else
    {
        for (uint16_t i = 0; i < NWORDS_FIELD - 1; i++)
            printf("%#018lx, ", ((digit_t *)el)[i]);
        printf("%#018lx", ((digit_t *)el)[NWORDS_FIELD - 1]);
    }
}

void print_f2elm(felm_t *el)
{
    // printf("{");
    print_felm_hex(&el[0],false);
    print_felm_hex(&el[1],true);
    // printf("},");
    printf("\n");
}

void load_f2elm(felm_t *target, unsigned long long *in, int len)
{
    for (int i = 0; i < len; i++)
        ((digit_t *)target)[(i / NWORDS_FIELD) * NWORDS_FIELD + (i % NWORDS_FIELD)] = in[i];
}

void load_E(instance_t *inst, CurveAndPointsSIDH *E0, CurveAndPointsSIDH *E1)
{
    inst->E[0] = *E0;
    inst->E[1] = *E1;
}

void print_E(instance_t *inst)
{
    print_f2elm(inst->E[0].a24);
    print_f2elm(inst->E[0].xp);
    print_f2elm(inst->E[0].xq);
    print_f2elm(inst->E[0].xpq);
    printf(" at %llx\n", (long long unsigned int)&(inst->E[0]));
    print_f2elm(inst->E[1].a24);
    print_f2elm(inst->E[1].xp);
    print_f2elm(inst->E[1].xq);
    print_f2elm(inst->E[1].xpq);
    printf(" at %llx\n", (long long unsigned int)&(inst->E[1]));
}

CurveAndPointsSIDH *digit_t_to_CurveAndPointsSIDH_ptr(digit_t *ptr)
{
    return (CurveAndPointsSIDH *)ptr;
}
#endif

