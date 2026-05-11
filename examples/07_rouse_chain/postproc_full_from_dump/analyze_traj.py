#!/usr/bin/env python3
import argparse
import numpy as np


def read_rouse_trajectory(filename, first_type=2, last_type=3):
    """
    Read a LAMMPS trajectory dump containing the full Rouse-chain system.

    Required dump columns:
        id mol type xu yu zu

    Returns
    -------
    steps : array, shape (nframes,)
    com   : array, shape (nframes, nchains, 3)
        Chain center-of-mass trajectory.
    rg2   : array, shape (nframes, nchains)
        Chain radius-of-gyration squared.
    ree   : array, shape (nframes, nchains, 3)
        Chain end-to-end vectors.

    Assumes equal bead masses.
    """

    steps = []
    com_frames = []
    rg2_frames = []
    ree_frames = []

    reference_mols = None

    with open(filename, "r") as f:
        while True:
            line = f.readline()

            if not line:
                break

            line = line.strip()

            if not line:
                continue

            if line != "ITEM: TIMESTEP":
                raise ValueError(f"Expected ITEM: TIMESTEP, got: {line}")

            step = int(f.readline().strip())

            line = f.readline().strip()
            if line != "ITEM: NUMBER OF ATOMS":
                raise ValueError(f"Expected ITEM: NUMBER OF ATOMS, got: {line}")

            n_atoms = int(f.readline().strip())

            line = f.readline().strip()
            if not line.startswith("ITEM: BOX BOUNDS"):
                raise ValueError(f"Expected ITEM: BOX BOUNDS, got: {line}")

            # Skip box bounds.
            for _ in range(3):
                f.readline()

            header = f.readline().strip().split()

            if len(header) < 3 or header[0] != "ITEM:" or header[1] != "ATOMS":
                raise ValueError("Expected ITEM: ATOMS header.")

            cols = header[2:]

            required = ["mol", "type", "xu", "yu", "zu"]
            for name in required:
                if name not in cols:
                    raise ValueError(
                        f"Trajectory must contain column '{name}'. "
                        f"Available columns: {cols}"
                    )

            imol = cols.index("mol")
            itype = cols.index("type")
            ix = cols.index("xu")
            iy = cols.index("yu")
            iz = cols.index("zu")

            sums = {}
            sums_sq = {}
            counts = {}
            first = {}
            last = {}

            for _ in range(n_atoms):
                parts = f.readline().split()

                mol = int(float(parts[imol]))
                typ = int(float(parts[itype]))

                pos = np.array(
                    [
                        float(parts[ix]),
                        float(parts[iy]),
                        float(parts[iz]),
                    ],
                    dtype=float,
                )

                if mol not in sums:
                    sums[mol] = pos.copy()
                    sums_sq[mol] = float(np.dot(pos, pos))
                    counts[mol] = 1
                else:
                    sums[mol] += pos
                    sums_sq[mol] += float(np.dot(pos, pos))
                    counts[mol] += 1

                if typ == first_type:
                    first[mol] = pos.copy()
                elif typ == last_type:
                    last[mol] = pos.copy()

            mols = sorted(sums.keys())

            if reference_mols is None:
                reference_mols = mols
            elif mols != reference_mols:
                raise ValueError("Molecule list changed between frames.")

            nchains = len(reference_mols)

            com = np.zeros((nchains, 3), dtype=float)
            rg2 = np.zeros(nchains, dtype=float)
            ree = np.zeros((nchains, 3), dtype=float)

            for i, mol in enumerate(reference_mols):
                if mol not in first or mol not in last:
                    raise ValueError(
                        f"Molecule {mol} does not have both endpoint types. "
                        f"Check first_type={first_type}, last_type={last_type}."
                    )

                n = counts[mol]

                com_i = sums[mol] / n

                # Rg^2 = <r^2> - |<r>|^2
                rg2_i = sums_sq[mol] / n - float(np.dot(com_i, com_i))

                # Avoid tiny negative values from roundoff.
                rg2_i = max(rg2_i, 0.0)

                com[i, :] = com_i
                rg2[i] = rg2_i
                ree[i, :] = last[mol] - first[mol]

            steps.append(step)
            com_frames.append(com)
            rg2_frames.append(rg2)
            ree_frames.append(ree)

    if len(steps) == 0:
        raise ValueError(f"No trajectory frames found in {filename}")

    return (
        np.asarray(steps),
        np.asarray(com_frames),
        np.asarray(rg2_frames),
        np.asarray(ree_frames),
    )


def calculate_multi_origin_msd(com, max_lag_fraction=0.5, origin_stride=1):
    """
    Multi-time-origin chain-COM MSD.

    MSD(t) = < |R_cm(t0+t) - R_cm(t0)|^2 >_{t0, chains}
    """

    nframes, nchains, _ = com.shape

    max_lag = int(max_lag_fraction * (nframes - 1))
    max_lag = max(1, max_lag)

    origins_all = np.arange(0, nframes, origin_stride, dtype=int)

    lags = []
    msd = []
    counts = []

    for lag in range(max_lag + 1):
        origins = origins_all[origins_all + lag < nframes]

        if len(origins) == 0:
            break

        dr = com[origins + lag, :, :] - com[origins, :, :]
        dr2 = np.sum(dr * dr, axis=2)

        lags.append(lag)
        msd.append(np.mean(dr2))
        counts.append(dr2.size)

    return np.asarray(lags), np.asarray(msd), np.asarray(counts)


def calculate_ree_autocorrelation(ree, max_lag_fraction=0.5, origin_stride=1):
    """
    Multi-time-origin end-to-end vector autocorrelation.

    Cee(t) =
        < Ree(t0+t) . Ree(t0) >_{t0, chains} / < Ree^2 >
    """

    nframes, nchains, _ = ree.shape

    max_lag = int(max_lag_fraction * (nframes - 1))
    max_lag = max(1, max_lag)

    origins_all = np.arange(0, nframes, origin_stride, dtype=int)

    lags = []
    c_raw = []
    counts = []

    for lag in range(max_lag + 1):
        origins = origins_all[origins_all + lag < nframes]

        if len(origins) == 0:
            break

        v0 = ree[origins, :, :]
        vt = ree[origins + lag, :, :]

        dots = np.sum(vt * v0, axis=2)

        lags.append(lag)
        c_raw.append(np.mean(dots))
        counts.append(dots.size)

    lags = np.asarray(lags)
    c_raw = np.asarray(c_raw)
    counts = np.asarray(counts)

    c_norm = c_raw / c_raw[0]

    return lags, c_raw, c_norm, counts


def fit_diffusion(time, msd, counts, tmin=None, tmax=None):
    """
    Fit:

        MSD(t) = 6 D t + constant

    using weighted linear regression.
    """

    if tmin is None:
        tmin = 0.10 * time[-1]

    if tmax is None:
        tmax = 0.70 * time[-1]

    mask = (time >= tmin) & (time <= tmax)

    if np.count_nonzero(mask) < 3:
        raise ValueError("Too few points in MSD fit window.")

    coeff = np.polyfit(
        time[mask],
        msd[mask],
        deg=1,
        w=np.sqrt(counts[mask]),
    )

    slope = coeff[0]
    intercept = coeff[1]

    D = slope / 6.0

    return D, slope, intercept, tmin, tmax, np.count_nonzero(mask)


def fit_longest_relaxation_time(time, c_norm, cmin=0.15, cmax=0.65):
    """
    Fit the long-time decay:

        Cee(t) ~ A exp(-t/tau_R)

    using the region cmin < Cee < cmax.
    """

    mask = (
        (time > 0.0)
        & np.isfinite(c_norm)
        & (c_norm > cmin)
        & (c_norm < cmax)
    )

    if np.count_nonzero(mask) < 3:
        return np.nan, np.nan, np.nan, 0

    slope, intercept = np.polyfit(time[mask], np.log(c_norm[mask]), deg=1)

    if slope >= 0.0:
        return np.nan, slope, intercept, np.count_nonzero(mask)

    tau = -1.0 / slope

    return tau, slope, intercept, np.count_nonzero(mask)


def main():
    parser = argparse.ArgumentParser(
        description="Analyze Brownian Rouse-chain trajectory from one LAMMPS dump file."
    )

    parser.add_argument("--traj", default="traj_rouse.lammpstrj")

    parser.add_argument("--dt", type=float, default=0.002)
    parser.add_argument("--nbeads", type=int, default=32)
    parser.add_argument("--temperature", type=float, default=1.0)
    parser.add_argument("--gamma", type=float, default=10.0)
    parser.add_argument("--b", type=float, default=1.0)

    parser.add_argument("--first-type", type=int, default=2)
    parser.add_argument("--last-type", type=int, default=3)

    parser.add_argument("--skip-fraction", type=float, default=0.2)

    parser.add_argument("--msd-max-lag-fraction", type=float, default=0.5)
    parser.add_argument("--msd-origin-stride", type=int, default=1)
    parser.add_argument("--msd-fit-tmin", type=float, default=None)
    parser.add_argument("--msd-fit-tmax", type=float, default=None)

    parser.add_argument("--corr-max-lag-fraction", type=float, default=0.5)
    parser.add_argument("--corr-origin-stride", type=int, default=1)
    parser.add_argument("--corr-fit-cmin", type=float, default=0.15)
    parser.add_argument("--corr-fit-cmax", type=float, default=0.65)

    args = parser.parse_args()

    # -------------------------
    # Read trajectory
    # -------------------------

    steps, com, rg2_chain, ree = read_rouse_trajectory(
        args.traj,
        first_type=args.first_type,
        last_type=args.last_type,
    )

    if len(steps) < 2:
        raise ValueError("Need at least two trajectory frames.")

    step_interval = int(np.median(np.diff(steps)))
    time_interval = step_interval * args.dt
    time_frames = steps * args.dt

    nframes = len(steps)
    nchains = com.shape[1]

    # -------------------------
    # Radius of gyration
    # -------------------------

    rg_chain = np.sqrt(rg2_chain)

    avg_rg_t = np.mean(rg_chain, axis=1)
    avg_rg2_t = np.mean(rg2_chain, axis=1)

    skip = int(args.skip_fraction * nframes)

    avg_rg = np.mean(avg_rg_t[skip:])
    avg_rg2 = np.mean(avg_rg2_t[skip:])
    sqrt_avg_rg2 = np.sqrt(avg_rg2)

    np.savetxt(
        "o.avg_rg.dat",
        np.column_stack([time_frames, avg_rg_t, avg_rg2_t]),
        header="time  <Rg>_chains  <Rg^2>_chains",
    )

    # -------------------------
    # End-to-end size
    # -------------------------

    ree2_t = np.mean(np.sum(ree * ree, axis=2), axis=1)
    avg_ree2 = np.mean(ree2_t[skip:])

    np.savetxt(
        "o.avg_ree2.dat",
        np.column_stack([time_frames, ree2_t]),
        header="time  <Ree^2>_chains",
    )

    # -------------------------
    # Multi-time-origin chain-COM MSD
    # -------------------------

    lags_msd, msd, msd_counts = calculate_multi_origin_msd(
        com,
        max_lag_fraction=args.msd_max_lag_fraction,
        origin_stride=args.msd_origin_stride,
    )

    time_msd = lags_msd * time_interval

    d_cm_fit, msd_slope, msd_intercept, fit_tmin, fit_tmax, nfit_msd = fit_diffusion(
        time_msd,
        msd,
        msd_counts,
        tmin=args.msd_fit_tmin,
        tmax=args.msd_fit_tmax,
    )

    np.savetxt(
        "o.avg_cm_msd_multi_origin.dat",
        np.column_stack([time_msd, msd, msd_counts]),
        header="time  <MSD_cm>_multi_origin  counts",
    )

    # -------------------------
    # End-to-end autocorrelation
    # -------------------------

    lags_corr, c_raw, c_norm, corr_counts = calculate_ree_autocorrelation(
        ree,
        max_lag_fraction=args.corr_max_lag_fraction,
        origin_stride=args.corr_origin_stride,
    )

    time_corr = lags_corr * time_interval

    tau_fit, corr_slope, corr_intercept, nfit_corr = fit_longest_relaxation_time(
        time_corr,
        c_norm,
        cmin=args.corr_fit_cmin,
        cmax=args.corr_fit_cmax,
    )

    np.savetxt(
        "o.ree_autocorr.dat",
        np.column_stack([time_corr, c_raw, c_norm, corr_counts]),
        header="time  <Ree(t).Ree(0)>  Cee_normalized  counts",
    )

    # -------------------------
    # Analytical references
    # -------------------------

    N = args.nbeads
    T = args.temperature
    gamma = args.gamma
    b = args.b

    d_cm_ref = T / (gamma * N)

    ree2_ref = (N - 1.0) * b**2

    rg2_ref = ((N**2 - 1.0) * b**2) / (6.0 * N)
    rg_ref = np.sqrt(rg2_ref)

    tau_rouse_continuum = gamma * b**2 * N**2 / (3.0 * np.pi**2 * T)

    tau_rouse_discrete = (
        gamma * b**2
        / (12.0 * T * np.sin(np.pi / (2.0 * N)) ** 2)
    )

    tau_from_rg_d = (2.0 / np.pi**2) * avg_rg2 / d_cm_fit

    # -------------------------
    # Print summary
    # -------------------------

    print("")
    print("Brownian Rouse-chain analysis from one trajectory")
    print("-------------------------------------------------")
    print(f"Input trajectory                         = {args.traj}")
    print(f"Number of frames                         = {nframes}")
    print(f"Number of chains                         = {nchains}")
    print(f"Number of beads per chain N              = {N}")
    print(f"Trajectory dump interval                 = {step_interval} timesteps")
    print(f"Trajectory dump interval                 = {time_interval:.6f} time units")
    print("")

    print("Model parameters")
    print("----------------")
    print(f"Temperature T                            = {T}")
    print(f"Friction gamma                           = {gamma}")
    print(f"Statistical segment length b             = {b}")
    print("")

    print("Radius of gyration")
    print("------------------")
    print(f"<Rg> averaged over chains and time        = {avg_rg:.6f}")
    print(f"<Rg^2> averaged over chains and time      = {avg_rg2:.6f}")
    print(f"sqrt(<Rg^2>)                              = {sqrt_avg_rg2:.6f}")
    print(f"Analytical <Rg^2>                         = {rg2_ref:.6f}")
    print(f"Analytical sqrt(<Rg^2>)                   = {rg_ref:.6f}")
    print("")

    print("End-to-end vector size")
    print("----------------------")
    print(f"<Ree^2> averaged over chains and time     = {avg_ree2:.6f}")
    print(f"Analytical <Ree^2>                        = {ree2_ref:.6f}")
    print("")

    print("Chain-COM diffusion: multi-time-origin")
    print("--------------------------------------")
    print(f"MSD fit window                           = {fit_tmin:.6f} to {fit_tmax:.6f}")
    print(f"Number of MSD fit points                 = {nfit_msd}")
    print(f"Fitted D_cm                              = {d_cm_fit:.6e}")
    print(f"Analytical D_cm                          = {d_cm_ref:.6e}")
    print("")

    print("End-to-end autocorrelation")
    print("--------------------------")
    print(f"Cee(0)                                    = {c_norm[0]:.6f}")
    print(f"Fitted tau_R from Cee(t)                  = {tau_fit:.6f}")
    print(f"Number of Cee fit points                  = {nfit_corr}")
    print(f"Fit window                                = {args.corr_fit_cmin} < Cee < {args.corr_fit_cmax}")
    print(f"Continuum Rouse tau_R                     = {tau_rouse_continuum:.6f}")
    print(f"Discrete-chain tau_R                      = {tau_rouse_discrete:.6f}")
    print(f"tau_R estimated from 2 Rg^2/(pi^2 D)      = {tau_from_rg_d:.6f}")
    print("")

    print("Output files written:")
    print("  avg_rg.dat")
    print("  avg_ree2.dat")
    print("  avg_cm_msd_multi_origin.dat")
    print("  ree_autocorr.dat")
    print("")


if __name__ == "__main__":
    main()
