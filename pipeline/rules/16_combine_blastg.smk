# rules/16_combine_blastg.smk

OUTPUT_DIR = config["paths"]["output"]

rule combine_blastg:
    """
    Combine ghost_data.csv and ghost_sum.csv from all chunk_* directories.
    """
    input:
        ghost_data = expand(f"{OUTPUT_DIR}/15_blastg/chunk_{{i:02d}}/ghost_data.csv", i=range(1, N_CHUNKS+1)),
        ghost_sum  = expand(f"{OUTPUT_DIR}/15_blastg/chunk_{{i:02d}}/ghost_sum.csv", i=range(1, N_CHUNKS+1))
    output:
        combined_data = f"{OUTPUT_DIR}/16_combine_blastg/ghost_data.csv",
        combined_sum  = f"{OUTPUT_DIR}/16_combine_blastg/ghost_sum.csv"
    conda:
        "envs/environment.yaml"
    params:
        combine_blastg_dir = f"{OUTPUT_DIR}/16_combine_blastg"
    run:
        import os
        import pandas as pd

        os.makedirs(params.combine_blastg_dir, exist_ok=True)

        # Combine ghost_data.csv
        dfs = [pd.read_csv(f) for f in input.ghost_data]
        pd.concat(dfs, ignore_index=True).to_csv(output.combined_data, index=False)

        # Combine ghost_sum.csv
        dfs = [pd.read_csv(f) for f in input.ghost_sum]
        pd.concat(dfs, ignore_index=True).to_csv(output.combined_sum, index=False)
