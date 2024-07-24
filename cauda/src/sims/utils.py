from micom import Community, data
from micom.qiime_formats import load_qiime_model_db, load_qiime_medium
import numpy as np
import pandas as pd
import warnings


def calc_kos(medium_series, metab_df, indices):
  """
  Calculate the effects of single member knockouts on a community.

  Parameters:
    medium_series (pandas Series) : A pandas series containing the base medium formulation for the simulation.
    metab_df (pandas DataFrame) : The metabolic database metadata as a pandas dataframe.
    indices (list of ints) : A list of indices from `metab_df` to use in the community.
  """
  tax = metab_df.loc[indices,]
  tax['id'] = tax['id'].apply(str)
    
  try:
    com = Community(tax,
                      progress=False)
    com.medium = medium_series
    ko = com.knockout_taxa(fraction=1.0, 
                             method="relative change", 
                             diag=False,
                             progress=False)
 
  # convert dataframe from wide to long format
    ko['knockout'] = ko.index
    ko = ko.melt(id_vars=["knockout"], 
               var_name="kept", 
               value_name="kept_relative_change")
    ko = ko.dropna()
    ko['error'] = 'None'

  # log errors for failed knockout
  except Exception as e:
    ko = pd.DataFrame({"knockout":list(tax.id),
                       "kept":list(tax.id)[::-1],
                       "kept_relative_change":[np.nan]*2})
    ko['error'] = e

  ko['SampleID'] = "_".join(sorted(list(tax.id)))

  return(ko)


def interp_fluxes(medium_1, medium_2, n, 
                       flux_col_1='flux', flux_col_2='flux'):
    """
    Interpolates fluxes between two chemical environments.
    
    Parameters:
    medium_1 (pd.DataFrame): DataFrame containing the fluxes for the medium_1.
    medium_2 (pd.DataFrame): DataFrame containing the fluxes for the medium_2.
    n (float): Interpolation factor (0 <= n <= 1).
    
    Returns:
    pd.Series: Interpolated fluxes.
    """
    if not (0 <= n <= 1):
        raise ValueError("n must be between 0 and 1")
    
    # Extract the flux series from the dataframes
    medium_1_series = medium_1[flux_col_1]
    medium_2_series = medium_2[flux_col_2]
    
    # Combine the indices of both series and remove duplicates
    all_rxns = list(set(medium_1_series.index).union(set(medium_2_series.index)))
    
    # Ensure indices are unique and create DataFrame with both series
    medium_1_series = medium_1_series[~medium_1_series.index.duplicated(keep='first')]
    medium_2_series = medium_2_series[~medium_2_series.index.duplicated(keep='first')]
    
    combined_df = pd.DataFrame(index=all_rxns)
    combined_df['medium_1_flux'] = medium_1_series.reindex(combined_df.index)
    combined_df['medium_2_flux'] = medium_2_series.reindex(combined_df.index)
    
    # Fill NaN values with 0 (assuming missing flux is 0)
    combined_df = combined_df.fillna(0)
    
    # Linearly interpolate the fluxes
    interpolated_flux = (1 - n) * combined_df['medium_1_flux'] + n * combined_df['medium_2_flux']
    return interpolated_flux