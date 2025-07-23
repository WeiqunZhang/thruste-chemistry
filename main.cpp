#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

#include "mechanism.H"

using namespace amrex;

namespace {
    constexpr int qrho = 0;
    constexpr int qtemp = 1;
    constexpr int qrhoy = 2;
}

void burn (MultiFab& q)
{
    auto const& ma = q.arrays();

    int nsteps = 10;
    Real dt = 1.e-12;
    for (int istep = 0; istep < nsteps; ++istep) {
        ParallelFor(q, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            auto const& a = ma[b].cellData(i,j,k);
            Real Yt[NUM_SPECIES];
            for (int n = 0; n < NUM_SPECIES; ++n) {
                Yt[n] = a[qrhoy+n] * (1._rt/a[qrho]);
            }

            Real wdot[NUM_SPECIES];
            CKWYR(a[qrho], a[qtemp], Yt, wdot);

            CKWT(Yt); // Yt now stores molecular weight;
            for (int n = 0; n < NUM_SPECIES; ++n) {
                a[qrhoy+n] += dt * wdot[n] * Yt[n];
            }
        });
        Gpu::streamSynchronize();
    }
}

void main_main ();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    main_main();
    amrex::Finalize();
}

void main_main ()
{
    int n_cell = 128;
    int max_grid_size = 64;
    Real problo = -0.8;
    Real probhi =  0.8;
    Real rfire  = 0.01;
    int nsteps;
    {
        ParmParse pp;
        pp.query("n_cell", n_cell);
        pp.query("max_grid_size", max_grid_size);
    }

    Box domain(IntVect(0),IntVect(n_cell-1));
    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    DistributionMapping dm{ba};

    MultiFab q(ba, dm, qrhoy+NUM_SPECIES, 0);

    Real dx = (probhi - problo) / Real(n_cell);
    auto const& ma = q.arrays();
    ParallelFor(q, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        auto x = problo + dx*(i+Real(0.5));
        auto y = problo + dx*(j+Real(0.5));
        auto z = problo + dx*(k+Real(0.5));
        auto r = std::sqrt(x*x + y*y + z*z);
        Real ru, ruc, patm;
        CKRP(ru, ruc, patm);

        auto Pt = patm;
        auto Tt = 300._rt;

        Real Xt[NUM_SPECIES];
        for (int n = 0; n < NUM_SPECIES; ++n) {
            Xt[n] = 1.e-10_rt;
        }
        Xt[H2_ID] = 0.10_rt;
        Xt[O2_ID] = 0.25_rt;

        auto expfac = std::exp(-Math::powi<2>(r/rfire));
        Pt += 0.1_rt*patm * expfac;
        Tt += 1100._rt * expfac;
        Xt[H2_ID] +=  0.025_rt * expfac;
        Xt[O2_ID] += -0.050_rt * expfac;
        Xt[N2_ID] = 1._rt - Xt[H2_ID] - Xt[O2_ID] - (NUM_SPECIES-3)*1.e-10_rt;

        Real Yt[NUM_SPECIES];
        CKXTY(Xt, Yt);
        Real rhot;
        CKRHOY(Pt,Tt,Yt,rhot);

        ma[b](i,j,k,qrho) = rhot;
        ma[b](i,j,k,qtemp) = Tt;
        for (int n = 0; n < NUM_SPECIES; ++n) {
            ma[b](i,j,k,qrhoy+n) = rhot*Yt[n];
        }
    });

    burn(q);
}
