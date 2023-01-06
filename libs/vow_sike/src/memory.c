#include <string.h>
#include "memory.h"
#include "types/triples.h"
#include "settings.h"
#include "triples.h"

int initialize_shared_memory(shared_state_t *S)
{
    if ((S->memory = calloc(S->MEMORY_SIZE, sizeof(trip_t))) == NULL) // triples
        return EXIT_FAILURE;
    for (uint64_t i = 0; i < S->MEMORY_SIZE; i++)
    {
        if ((S->memory[i].current_state.words = calloc(S->NWORDS_STATE, sizeof(digit_t))) == NULL) /* Only words since > bytes */
            return EXIT_FAILURE;
        if ((S->memory[i].initial_state.words = calloc(S->NWORDS_STATE, sizeof(digit_t))) == NULL)
            return EXIT_FAILURE;
        S->memory[i].current_steps = 0;
    }
    return EXIT_SUCCESS;
}

void initialize_private_memory(shared_state_t *S, private_state_t *private_state)
{
    (void)S;
    (void)private_state;
}

void cleanup_shared_memory(shared_state_t *shared_state)
{
    for (unsigned int i = 0; i < shared_state->MEMORY_SIZE; i++)
    {
        free_trip(&shared_state->memory[i]);
    }
    free(shared_state->memory);
}

void cleanup_private_memory(private_state_t *private_state)
{
    // Do nothing
    (void)private_state;
}

/*
 * Reads triple from memory at specified address.
 * If no triple is found there, returns false, else true.
 */
extern inline bool read_from_memory(trip_t **t, shared_state_t *S, private_state_t *private_state, uint64_t address)
{
    (void)private_state;
    *t = &S->memory[address];
    return (*t)->current_steps > 0; // if no steps, then memory location was empty
}

extern inline void write_to_memory(trip_t *t, shared_state_t *S, private_state_t *private_state, uint64_t address)
{
    (void)private_state;
    copy_trip(&S->memory[address], t, private_state->NWORDS_STATE);
}

void fix_overflow(st_t *s, const uint64_t nbytes_state, const uint64_t nbits_overflow)
{
    if (nbits_overflow != 0)
        s->bytes[nbytes_state - 1] &= (0xFF >> (8 - nbits_overflow));
}