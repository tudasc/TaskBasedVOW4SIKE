#pragma once
#include "types/state.h"

void nobiggie_resync_do_resync(shared_state_t *S, private_state_t *private_state);
void print_all_threads(shared_state_t *S, private_state_t *private_state, unsigned char *buffer);
void stakhanovist_resync_do_resync(shared_state_t *S, private_state_t *private_state);
bool stakhanovist_resync_should_resync(shared_state_t *S, private_state_t *private_state);
void windowed_resync(shared_state_t *S, private_state_t *private_state);
void windowed_resync_do_resync(shared_state_t *S, private_state_t *private_state);
bool windowed_resync_should_resync(shared_state_t *S, private_state_t *private_state);