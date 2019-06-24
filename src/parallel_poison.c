#include "mpi.h"

#include "parallel_poisson.h"
#include "poissonNonlinPeriodic.h"
#include "poissonNonlinNonPeriodic.h"
#include "zBar.h"

void get_uniform_Te(double **Te, int *Nx, double *T0){
    for(int i = 0; i < *Nx; i++){
        (*Te)[i] = *T0;
    }
}

void get_bc_Te(double **Te_arr_allranks, int *Nx, int *order, double **T0_bcs){
    for(int i = 0; i < *order; i++){
        (*Te_arr_allranks)[i] = (*T0_bcs)[0];
        (*Te_arr_allranks)[(*Nx) + i] = (*T0_bcs)[1];
    }
}

void get_ramp_Te(double **Te, int *Nx, double *Te_start, double *Te_end, double *t, double *tfinal){
    double rampTemp;
    if(*t > *tfinal){
        rampTemp = (*Te_end);
    }else{
        rampTemp = (*Te_start) + (*t)*((*Te_end) - (*Te_start))/(*tfinal);
    }
    for(int l = 0; l < *Nx; l++){
        (*Te)[l] = rampTemp;
    }
}

void get_ramp_bc_Te(double **Te_arr_allranks, int *Nx, int *order, double *T_start, 
                    double *T_end, double *t, double *tfinal)
{
    double rampTemp;
    if(*t > *tfinal){
        rampTemp = (*T_end);
    }else{
        rampTemp = (*T_start) + (*t)*((*T_end)-(*T_start))/(*tfinal);
    }
    for(int l = 0; l < *order; l++){
        (*Te_arr_allranks)[l] = rampTemp;
        (*Te_arr_allranks)[(*Nx) + l] = rampTemp;
    }
}


void periodic_poisson_solver(MPI_Comm comm, int *rank, int *numRanks, MPI_Status *status, int *Nx, int **Nx_ranks, int *Nx_rank,
                    int *ecouple, int *poissFlavor, int *ionFix, int *nspec, double **Z_max, double *Te_ref, 
                    double *Te_start, int *order, double *dx, double *Lx, double *t, double *tfinal, double **Te_arr,
                    double **Te_arr_allranks, double **T0_oned, double ***n_oned, double ***Z_oned, 
                    double **source, double **source_buf, double **source_allranks, double **PoisPot, double **PoisPot_allranks)
{
    //set electron temp. This assumes that dynamic electrons have the local mixture temp,
    //i.e. equilibrate extremely rapidly.
    if(*ecouple == 1){
        //electrons only in background
        if(*Te_start != *Te_ref){
            get_ramp_Te(Te_arr, Nx_rank, Te_start, Te_ref, t, tfinal);
        }else{
            get_uniform_Te(Te_arr, Nx_rank, Te_ref);
        }
    }else{
        //local equilibrium dynamic electrons
        *Te_arr = *T0_oned;
    }
    //Calculate ionization data, store source terms in buffer.
    for(int l = 0; l < *Nx_rank; l++){
        if((*ecouple != 2) && (*ionFix != 1)){
            zBarFunc2(*nspec, (*Te_arr)[l], *Z_max, (*n_oned)[l], (*Z_oned)[l]);
        }else{
            (*Z_oned)[l] = (*Z_max);
        }
        (*source)[l] = 0.0;
        if(*ecouple != 2){
            (*source)[l] += (*Z_oned)[l][0]*(*n_oned)[l][0];
        }else{
            (*source)[l] -= (*Z_oned)[l][0]*(*n_oned)[l][0];
        }
        for(int i = 1; i < *nspec; i++){
            (*source)[l] += (*Z_oned)[l][i];
        }
    }
    //Set up source/RHS array for Poisson solve
    if(*rank == 0){
        for(int l = 0; l < *Nx_rank; l++){
            (*source_allranks)[l] = (*source)[l];
            (*Te_arr_allranks)[l] = (*Te_arr)[l];
        }

        if(*numRanks > 1){
            int rankOffset = *Nx_rank;
            //get source value and electron Temp from all other ranks.
            for(int rankCounter = 1; rankCounter < *numRanks; rankCounter++){
                MPI_Recv(*source_buf, 2 * (*Nx_ranks)[rankCounter], MPI_DOUBLE,
                         rankCounter, 200 + rankCounter, comm, status);
                for(int l = 0; l < (*Nx_ranks)[rankCounter]; l++){
                    (*source_allranks)[l + rankOffset] = (*source_buf)[0 + 2*l];
                    (*Te_arr_allranks)[l + rankOffset] = (*source_buf)[1 + 2*l];
                }
                rankOffset += (*Nx_ranks)[rankCounter];
            }
        }
    }else{
        //send rhs value to rank 0.
        for(int l = 0; l < *Nx_rank; l++){
            (*source_buf)[0 + 2*l] = (*source)[l];
            (*source_buf)[1 + 2*l] = (*Te_arr)[l];
        }
        MPI_Send(*source_buf, 2 * (*Nx_rank), MPI_DOUBLE, 0, 200 + (*rank), comm);
    }

    //Wait until all source stuff is sent.
    MPI_Barrier(comm);
    
    //Rank 0 performs the solve.
    if(*rank == 0){
        if(*poissFlavor == 0){
            //No E-field
            for(int l = 0; l < *Nx; l++){
                (*PoisPot_allranks)[l] = 0.0;
            }
        }else if(*poissFlavor == 11){
            //Linear Yukawa
            PoissLinPeriodic1D(*Nx, *source_allranks, *dx, *Lx, *PoisPot_allranks, *Te_arr_allranks);
        }else if(*poissFlavor == 12){
            //Nonlinear Yukawa
            PoissNonlinPeriodic1D(*Nx, *source_allranks, *dx, *Lx, *PoisPot_allranks, *Te_arr_allranks);
        }else if(*poissFlavor == 21){
            //Linear Thomas-Fermi
            PoissLinPeriodic1D_TF(*Nx, *source_allranks, *dx, *Lx, *PoisPot_allranks, *Te_arr_allranks);
        }else if(*poissFlavor == 22){
            PoissNonlinPeriodic1D_TF(*Nx, *source_allranks, *dx, *Lx, *PoisPot_allranks, *Te_arr_allranks);
        }
    }

    //Distribute back to the other ranks.
    if(*rank == 0){
        int rankOffset = 0;

        //Set local ghost cells.
        if(*order == 1){
            (*PoisPot)[0] = (*PoisPot_allranks)[*Nx -1];
            if(*numRanks == 1){
                (*PoisPot)[(*Nx_rank) + 1] = (*PoisPot_allranks)[0];
            }else{
                (*PoisPot)[(*Nx_rank) + 1] = (*PoisPot_allranks)[(*Nx_rank)];
            }
        }else if(*order == 2){
            (*PoisPot)[0] = (*PoisPot_allranks)[(*Nx)-2];
            (*PoisPot)[1] = (*PoisPot_allranks)[(*Nx)-1];
            if(*numRanks == 1){
                (*PoisPot)[(*Nx_rank) + 2] = (*PoisPot_allranks)[0];
                (*PoisPot)[(*Nx_rank) + 3] = (*PoisPot_allranks)[1];
            }else{
                (*PoisPot)[(*Nx_rank) + 2] = (*PoisPot_allranks)[(*Nx_rank)];
                (*PoisPot)[(*Nx_rank) + 3] = (*PoisPot_allranks)[(*Nx_rank) + 1];
            }
        }

        //Set main body of PoisPot
        for(int l = 0; l < *Nx_rank; l++){
            (*PoisPot)[l + (*order)] = (*PoisPot_allranks)[l];
        }

        if(*numRanks > 1){
            rankOffset = *Nx_rank;
            for(int rankCounter = 1; rankCounter < (*numRanks)-1; rankCounter++){
                if (*order == 1) {
                    (*source_buf)[0] = (*PoisPot_allranks)[rankOffset - 1];
                    (*source_buf)[(*Nx_ranks)[rankCounter] + 1] = 
                    (*PoisPot_allranks)[rankOffset + (*Nx_ranks)[rankCounter]];
                 } else if (*order == 2) {
                     (*source_buf)[0] = (*PoisPot_allranks)[rankOffset - 2];
                     (*source_buf)[1] = (*PoisPot_allranks)[rankOffset - 1];
                     (*source_buf)[(*Nx_ranks)[rankCounter] + 2] =
                     (*PoisPot_allranks)[rankOffset + (*Nx_ranks)[rankCounter]];
                     (*source_buf)[(*Nx_ranks)[rankCounter] + 3] =
                     (*PoisPot_allranks)[rankOffset + (*Nx_ranks)[rankCounter] + 1];
                }

                //Set main body of PoisPot on rankCounter
                for(int l = 0; l < (*Nx_ranks)[rankCounter]; l++){
                    (*source_buf)[l + (*order)] = (*PoisPot_allranks)[rankOffset + l];
                }

                MPI_Send(*source_buf, (*Nx_ranks)[rankCounter] + 2*(*order), MPI_DOUBLE,
                         rankCounter, rankCounter, comm);
                rankOffset += (*Nx_ranks)[rankCounter];
            }

            //Deal with periodic BC for rightmost rank

            //set local ghost cells
            if(*order == 1){
                (*source_buf)[0] = (*PoisPot_allranks)[rankOffset - 1];
                (*source_buf)[(*Nx_ranks)[(*numRanks) - 1] + 1] = (*PoisPot_allranks)[0];
            } else if(*order == 2) {
                (*source_buf)[0] = (*PoisPot_allranks)[rankOffset - 2];
                (*source_buf)[1] = (*PoisPot_allranks)[rankOffset - 1];
                (*source_buf)[(*Nx_ranks)[(*numRanks) - 1] + 2] = (*PoisPot_allranks)[0];
                (*source_buf)[(*Nx_ranks)[(*numRanks) - 1] + 3] = (*PoisPot_allranks)[1];
            }

            for(int l = 0; l < (*Nx_ranks)[(*numRanks)-1]; l++){
                (*source_buf)[l + (*order)] = (*PoisPot_allranks)[rankOffset + l];
            }
            
            MPI_Send(*source_buf, (*Nx_ranks)[(*numRanks)-1] + 2*(*order), MPI_DOUBLE,
                     (*numRanks)-1, (*numRanks)-1, comm);
        }
    }else{
        //Get potential from rank 0
        MPI_Recv(*source_buf, *Nx_rank + 2*(*order), MPI_DOUBLE, 0,  *rank, comm, status);

        //set Poispot
        for(int l = 0; l < (*Nx_rank) + 2*(*order); l++){
            (*PoisPot)[l] = (*source_buf)[l];
        }
    }
}

void nonperiodic_poisson_solver(MPI_Comm comm, int *rank, int *numRanks, MPI_Status *status, int *Nx, int **Nx_ranks, int *Nx_rank,
                    int *ecouple, int *poissFlavor, int *ionFix, int *nspec, double **Z_max, double *Te_ref, 
                    double *Te_start, int *order, double *dx, double *Lx, double *t, double *tfinal, double **Te_arr,
                    double **Te_arr_allranks, double **T0_oned, double ***n_oned, double ***Z_oned, 
                    double **source, double **source_buf, double **source_allranks, double **PoisPot, double **PoisPot_allranks,
                    double **T0_bcs, double ***n_bcs, double ***Z_bcs, double **Te_bcs)
{
    //set electron temp. This assumes that dynamic electrons have the local mixture temp,
    //i.e. equilibrate extremely rapidly.
    if(*ecouple == 1){
        //electrons only in background
        if(*Te_start != *Te_ref){
            get_ramp_Te(Te_arr, Nx_rank, Te_start, Te_ref, t, tfinal);
            if(*rank == 0){
                get_ramp_bc_Te(Te_arr_allranks, Nx, order, Te_start, Te_ref, t, tfinal);
            }
        }else{
            get_uniform_Te(Te_arr, Nx_rank, Te_ref);
            if(*rank == 0){
                for(int i = 0; i < *order; i++){
                    (*Te_arr_allranks)[i] = (*Te_ref);
                    (*Te_arr_allranks)[(*Nx) + i] = (*T0_bcs)[1];
                }
            }
        }
    }else{
        //local equilibrium dynamic electrons
        *Te_arr = *T0_oned;
        if(*rank == 0){
            get_bc_Te(Te_arr_allranks, Nx, order, T0_bcs);
        }
    }
    (*Te_bcs)[0] = (*Te_arr_allranks)[0];
    (*Te_bcs)[1] = (*Te_arr_allranks)[(*Nx)-1];
    //Calculate ionization data, store source terms in buffer.
    for(int l = 0; l < *Nx_rank; l++){
        if((*ecouple != 2) && (*ionFix != 1)){
            zBarFunc2(*nspec, (*Te_arr)[l], *Z_max, (*n_oned)[l], (*Z_oned)[l]);
        }else{
            (*Z_oned)[l] = (*Z_max);
        }
        (*source)[l] = 0.0;
        if(*ecouple != 2){
            (*source)[l] += (*Z_oned)[l][0]*(*n_oned)[l][0];
        }else{
            (*source)[l] -= (*Z_oned)[l][0]*(*n_oned)[l][0];
        }
        for(int i = 1; i < *nspec; i++){
            (*source)[l] += (*Z_oned)[l][i];
        }
    }
    //calculate ionization data for bcs, store in buffer.
    if(*rank == 0){
        if((*ecouple != 2) && (*ionFix != 1)){
            zBarFunc2(*nspec, (*Te_bcs)[0], *Z_max, (*n_bcs)[0], (*Z_bcs)[0]);
            zBarFunc2(*nspec, (*Te_bcs)[0], *Z_max, (*n_bcs)[1], (*Z_bcs)[1]);
        }else{
            (*Z_bcs)[0] = (*Z_max);
            (*Z_bcs)[1] = (*Z_max);
        }
        for(int l = 0; l < *order; l++){
            (*source_allranks)[l] = 0.0;
            (*source_allranks)[(*Nx) + l] = 0.0;
            if(*ecouple != 2){
                (*source_allranks)[l] += (*Z_bcs)[0][0]*(*n_bcs)[0][0];
                (*source_allranks)[(*Nx) + l] += (*Z_bcs)[1][0]*(*n_bcs)[1][0];
            }else{
                (*source_allranks)[l] -= (*Z_bcs)[0][0]*(*n_bcs)[0][0];
                (*source_allranks)[(*Nx) + l] -= (*Z_bcs)[1][0]*(*n_bcs)[1][0];
 
            }
            for(int i = 1; i < *nspec; i++){
                (*source_allranks)[l] += (*Z_bcs)[0][i]*(*n_bcs)[0][i];
                (*source_allranks)[(*Nx) + (*order) + l] +=  (*Z_bcs)[1][i]*(*n_bcs)[1][i];
            }
        }
    }
    //Set up source/RHS array for Poisson solve
    if(*rank == 0){
        for(int l = 0; l < *Nx_rank; l++){
            (*source_allranks)[l + (*order)] = (*source)[l];
            (*Te_arr_allranks)[l + (*order)] = (*Te_arr)[l];
        }

        if(*numRanks > 1){
            int rankOffset = *Nx_rank + *order;

            //get source value and electron Temp from all other ranks.
            for(int rankCounter = 1; rankCounter < *numRanks; rankCounter++){
                MPI_Recv(*source_buf, 2 * (*Nx_ranks)[rankCounter], MPI_DOUBLE,
                         rankCounter, 200 + rankCounter, comm, status);
                for(int l = 0; l < (*Nx_ranks)[rankCounter]; l++){
                    (*source_allranks)[l + rankOffset] = (*source_buf)[0 + 2*l];
                    (*Te_arr_allranks)[l + rankOffset] = (*source_buf)[1 + 2*l];
                }
                rankOffset += (*Nx_ranks)[rankCounter];
            }
        }
    }else{
        //send rhs value to rank 0.
        for(int l = 0; l < *Nx_rank; l++){
            (*source_buf)[0 + 2*l] = (*source)[l];
            (*source_buf)[1 + 2*l] = (*Te_arr)[l];
        }
        MPI_Send(*source_buf, 2 * (*Nx_rank), MPI_DOUBLE, 0, 200 + (*rank), comm);
    }

    //Wait until all source stuff is sent.
    MPI_Barrier(comm);
    
    //Rank 0 performs the solve.
    if(*rank == 0){
        if(*poissFlavor == 0){
            //No E-field
            for(int l = 0; l < *Nx; l++){
                (*PoisPot_allranks)[l] = 0.0;
            }
        }else if(*poissFlavor == 11){
            //Linear Yukawa
            PoissLinNonPeriodic1D(*Nx, order, *source_allranks, *dx, *Lx, *PoisPot_allranks, *Te_arr_allranks);
        }else if(*poissFlavor == 12){
            //Nonlinear Yukawa
            PoissNonlinNonPeriodic1D(*Nx, order, *source_allranks, *dx, *Lx, *PoisPot_allranks, *Te_arr_allranks);
        }else if(*poissFlavor == 21){
            //Linear Thomas-Fermi
            PoissLinNonPeriodic1D_TF(*Nx, order, *source_allranks, *dx, *Lx, *PoisPot_allranks, *Te_arr_allranks);
        }else if(*poissFlavor == 22){
            PoissNonlinNonPeriodic1D_TF(*Nx, order, *source_allranks, *dx, *Lx, *PoisPot_allranks, *Te_arr_allranks);
        }
    }

    //Distribute back to the other ranks.
    if(*rank == 0){
        int rankOffset = 0;

        //Set local ghost cells.
        if(*order == 1){
            (*PoisPot)[0] = (*PoisPot_allranks)[0];
            if(*numRanks == 1){
                (*PoisPot)[(*Nx_rank) + 1] = (*PoisPot_allranks)[(*Nx_rank) + 1];
            }else{
                (*PoisPot)[(*Nx_rank) + 1] = (*PoisPot_allranks)[(*Nx_rank) + 1];
            }
        }else if(*order == 2){
            (*PoisPot)[0] = (*PoisPot_allranks)[0];
            (*PoisPot)[1] = (*PoisPot_allranks)[1];
            if(*numRanks == 1){
                (*PoisPot)[(*Nx_rank) + 2] = (*PoisPot_allranks)[(*Nx_rank) + 2];
                (*PoisPot)[(*Nx_rank) + 3] = (*PoisPot_allranks)[(*Nx_rank) + 3];
            }else{
                (*PoisPot)[(*Nx_rank) + 2] = (*PoisPot_allranks)[(*Nx_rank) + 2];
                (*PoisPot)[(*Nx_rank) + 3] = (*PoisPot_allranks)[(*Nx_rank) + 3];
            }
        }

        //Set main body of PoisPot
        for(int l = (*order); l < *Nx_rank; l++){
            (*PoisPot)[l] = (*PoisPot_allranks)[l];
        }

        if(*numRanks > 1){
            rankOffset = *Nx_rank + (*order);
            for(int rankCounter = 1; rankCounter < (*numRanks)-1; rankCounter++){
                if (*order == 1) {
                    (*source_buf)[0] = (*PoisPot_allranks)[rankOffset - 1];
                    (*source_buf)[(*Nx_ranks)[rankCounter] + 1] = 
                    (*PoisPot_allranks)[rankOffset + (*Nx_ranks)[rankCounter]];
                 } else if (*order == 2) {
                     (*source_buf)[0] = (*PoisPot_allranks)[rankOffset - 2];
                     (*source_buf)[1] = (*PoisPot_allranks)[rankOffset - 1];
                     (*source_buf)[(*Nx_ranks)[rankCounter] + 2] =
                     (*PoisPot_allranks)[rankOffset + (*Nx_ranks)[rankCounter]];
                     (*source_buf)[(*Nx_ranks)[rankCounter] + 3] =
                     (*PoisPot_allranks)[rankOffset + (*Nx_ranks)[rankCounter] + 1];
                }

                //Set main body of PoisPot on rankCounter
                for(int l = 0; l < (*Nx_ranks)[rankCounter]; l++){
                    (*source_buf)[l + (*order)] = (*PoisPot_allranks)[rankOffset + l];
                }

                MPI_Send(*source_buf, (*Nx_ranks)[rankCounter] + 2*(*order), MPI_DOUBLE,
                         rankCounter, rankCounter, comm);
                rankOffset += (*Nx_ranks)[rankCounter];
            }

            //Deal with periodic BC for rightmost rank

            //set local ghost cells
            if(*order == 1){
                (*source_buf)[0] = (*PoisPot_allranks)[rankOffset - 1];
                (*source_buf)[(*Nx_ranks)[(*numRanks) - 1] + 1] = (*PoisPot_allranks)[(*Nx) + 2*(*order)-1];
            } else if(*order == 2) {
                (*source_buf)[0] = (*PoisPot_allranks)[rankOffset - 2];
                (*source_buf)[1] = (*PoisPot_allranks)[rankOffset - 1];
                (*source_buf)[(*Nx_ranks)[(*numRanks) - 1] + 2] = (*PoisPot_allranks)[(*Nx) + 2];
                (*source_buf)[(*Nx_ranks)[(*numRanks) - 1] + 3] = (*PoisPot_allranks)[(*Nx) + 3];
            }

            for(int l = 0; l < (*Nx_ranks)[(*numRanks)-1]; l++){
                (*source_buf)[l + (*order)] = (*PoisPot_allranks)[rankOffset + l];
            }
            
            MPI_Send(*source_buf, (*Nx_ranks)[(*numRanks)-1] + 2*(*order), MPI_DOUBLE,
                     (*numRanks)-1, (*numRanks)-1, comm);
        }
    }else{
        //Get potential from rank 0
        MPI_Recv(*source_buf, *Nx_rank + 2*(*order), MPI_DOUBLE, 0,  *rank, comm, status);

        //set Poispot
        for(int l = 0; l < (*Nx_rank) + 2*(*order); l++){
            (*PoisPot)[l] = (*source_buf)[l];
        }
    }
}



void poisson_solver(MPI_Comm comm, int *rank, int *numRanks, MPI_Status *status, int *Nx, int **Nx_ranks, int *Nx_rank,
                    int *ecouple, int *bcs, int *poissFlavor, int *ionFix, int *nspec, double **Z_max, double *Te_ref, 
                    double *Te_start, int *order, double *dx, double *Lx, double *t, double *tfinal, double **Te_arr,
                    double **Te_arr_allranks, double **T0_oned, double ***n_oned, double ***Z_oned, 
                    double **source, double **source_buf, double **source_allranks, double **PoisPot, double **PoisPot_allranks,
                    double **T0_bcs, double ***n_bcs, double ***Z_bcs, double **Te_bcs)
{
    if(*bcs == 0){
        periodic_poisson_solver(comm, rank, numRanks, status, Nx, Nx_ranks, Nx_rank,
                                ecouple, poissFlavor, ionFix, nspec, Z_max, Te_ref, 
                                Te_start, order, dx, Lx, t, tfinal, Te_arr, Te_arr_allranks,
                                T0_oned, n_oned, Z_oned, source, source_buf, source_allranks,
                                PoisPot, PoisPot_allranks);
    }else{
        nonperiodic_poisson_solver(comm, rank, numRanks, status, Nx, Nx_ranks, Nx_rank,
                    ecouple, poissFlavor, ionFix, nspec, Z_max, Te_ref, 
                    Te_start, order, dx, Lx, t, tfinal, Te_arr,
                    Te_arr_allranks, T0_oned, n_oned, Z_oned, 
                    source, source_buf, source_allranks, PoisPot, PoisPot_allranks,
                    T0_bcs, n_bcs, Z_bcs, Te_bcs);
    }
}