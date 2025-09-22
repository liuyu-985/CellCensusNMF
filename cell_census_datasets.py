# H matrix and cells metadata
import os
import pandas as pd
import cellxgene_census

output_embedding_file = "/mnt/projects/debruinz_project/yu_ting/adata_obsm.csv"
output_obs_file       = "/mnt/home/liuyu/adata_obs.csv"

os.makedirs(os.path.dirname(output_embedding_file), exist_ok=True)
os.makedirs(os.path.dirname(output_obs_file), exist_ok=True)

emb_header_written = False
obs_header_written = False

print("Opening Census (2023-12-15)...")
with cellxgene_census.open_soma(census_version="2023-12-15") as census:
    print("Retrieving Homo sapiens obs (tissue)...")
    obs_df = cellxgene_census.get_obs(census, "homo_sapiens", column_names=["tissue"])

    tissues = pd.Series(obs_df["tissue"]).dropna().unique().tolist()
    print(f"Found {len(tissues)} tissues to process.")

    for tissue in tissues:
        tissue_filter = tissue.replace("'", "\\'")
        print(f"\nProcessing tissue: {tissue}")

        try:
            adata = cellxgene_census.get_anndata(
                census,
                organism="homo_sapiens",
                measurement_name="RNA",
                obs_value_filter=f"tissue == '{tissue_filter}' and is_primary_data == True",
                obs_embeddings=["nmf"], 
            )

            
            adata.obs.to_csv(
                output_obs_file,
                mode="a",
                header=not obs_header_written,
                index=True,
            )
            obs_header_written = True
            print(f"Wrote obs rows for {tissue} → {output_obs_file}")

           
            emb = adata.obsm.get("nmf")
            if emb is not None:
                emb_df = pd.DataFrame(emb, index=adata.obs.index)

                
                if emb_df.columns.dtype.kind in ("i", "u"):  
                    emb_df.columns = [f"NMF_{i+1}" for i in range(emb_df.shape[1])]

                emb_df.to_csv(
                    output_embedding_file,
                    mode="a",
                    header=not emb_header_written,
                    index=True,
                )
                emb_header_written = True
                print(f"Wrote NMF embeddings for {tissue} → {output_embedding_file}")
            else:
                print(f"No 'nmf' embeddings available for {tissue}; skipped embeddings write.")

            
            del adata

        except Exception as e:
            print(f"Error processing tissue {tissue}: {e}")



# W matrix and genes metadata
import cellxgene_census
import pandas as pd
import os


census = cellxgene_census.open_soma(census_version="2023-12-15")


out_var_meta = "/mnt/projects/debruinz_project/yu_ting/adata_var_metadata.csv"
out_var_emb  = "/mnt/projects/debruinz_project/yu_ting/adata_var_embeddings.csv"
os.makedirs(os.path.dirname(out_var_meta), exist_ok=True)
os.makedirs(os.path.dirname(out_var_emb),  exist_ok=True)


meta_written = False
emb_written  = False


obs_df         = cellxgene_census.get_obs(census, "homo_sapiens", column_names=["tissue"])
unique_tissues = obs_df["tissue"].unique()


for tissue in unique_tissues:
    try:
        print(f"Processing tissue: {tissue}")
        adata = cellxgene_census.get_anndata(
            census,
            organism         = "homo_sapiens",
            measurement_name = "RNA",
            obs_value_filter = f"tissue == '{tissue}' and is_primary_data == True",
            var_embeddings   = ["nmf"]
        )

        
        var_meta = adata.var.copy()
        var_meta["tissue"] = tissue
        var_meta.to_csv(
            out_var_meta,
            mode   = "a",
            header = not meta_written,
            index  = True
        )
        meta_written = True

        
        if "nmf" in adata.varm:
            n_factors = adata.varm["nmf"].shape[1]
            var_emb = pd.DataFrame(
                adata.varm["nmf"],
                index   = adata.var.index,
                columns = [f"Factor{i+1}" for i in range(n_factors)]
            )
            var_emb["tissue"] = tissue
            var_emb.to_csv(
                out_var_emb,
                mode   = "a",
                header = not emb_written,
                index  = True
            )
            emb_written = True

    except Exception as e:
        print(f"Error for tissue {tissue}: {e}")

print(" Variable metadata and embeddings extraction complete.")



df = pd.read_csv(
    "/mnt/projects/debruinz_project/yu_ting/adata_var_embeddings.csv",
    index_col=0             
)
df1 = pd.read_csv(
    "/mnt/projects/debruinz_project/yu_ting/adata_var_metadata.csv",
    index_col=0             
)

df1_new = (
    df1
    .drop_duplicates(subset="feature_name", keep="first")
    .reset_index(drop=False)      
    .rename(columns={"index": "orig_row"})
)


rows_to_keep = df1_new["orig_row"].values
df_new   = df.iloc[rows_to_keep].reset_index(drop=True)


df1_new = df1_new.drop(columns=["orig_row", "tissue"], errors="ignore").reset_index(drop=True)
df_new  = df_new .drop(columns=["tissue"],          errors="ignore").reset_index(drop=True)


df1_new.to_csv(
    "/mnt/projects/debruinz_project/yu_ting/adata_var_metadata_unique.csv",
    index=True
)
df_new.to_csv(
    "/mnt/projects/debruinz_project/yu_ting/adata_var_embeddings_unique.csv",
    index=True
)

# transfer learning new data--transfer_learning_data.h5ad
https://cellxgene.cziscience.com/a6657cae-5daa-45cd-b1b4-cf08a07d3a7e.h5ad

