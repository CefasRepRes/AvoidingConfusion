#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np, plotly.graph_objects as go, plotly.offline as pyo
from bayesian_spence import (
    derive_dirichlet_posterior_parameters,
    extract_validation_error_model,
    load_validation_uncertainty_dataframe,
)
from ml_prediction_results_utilities import *
from plotting_predictions import plot_timeseries
def confusion(source,path):
 m=extract_validation_error_model(load_json_payload(source)); cs=m['class_order']; a=np.asarray(derive_dirichlet_posterior_parameters(m,cs)); z=a/a.sum(1,keepdims=True); fig=go.Figure(go.Heatmap(z=z,x=cs,y=cs,text=np.vectorize(lambda x:f'{x:.3f}')(z),texttemplate='%{text}',colorscale='Blues')); fig.update_layout(title='Classifier confusion matrix: posterior mean P(predicted | true)',xaxis_title='Predicted class',yaxis_title='True class'); pyo.plot(fig,filename=str(path),auto_open=False,include_plotlyjs='cdn')
def main():
 p=argparse.ArgumentParser(); p.add_argument('--input-json',required=True); p.add_argument('--validation-json',required=True); p.add_argument('--outdir',default='summary_timeseries_out'); p.add_argument('--target-class',default='fish_larvae'); p.add_argument('--mc-samples',type=int,default=250); p.add_argument('--mc-seed',type=int,default=42); p.add_argument('--stacked',action='store_true'); p.add_argument('--log-y',action='store_true'); p.add_argument('--visible-classes'); p.add_argument('--no-diagnostics',action='store_true'); a=p.parse_args(); o=Path(a.outdir); o.mkdir(parents=True,exist_ok=True); d,_=build_dataframe(load_inference_rows_from_source(a.input_json)); d.to_csv(o/'class_counts_timeseries.csv',index=False); write_long_counts_csv(d,o/'class_counts_timeseries_long.csv'); confusion(a.validation_json,o/'confusion_matrix_posterior_mean.html'); target=o/safe_part(a.target_class); target.mkdir(exist_ok=True); u=load_validation_uncertainty_dataframe(d,a.validation_json,mc_samples=a.mc_samples,mc_seed=a.mc_seed,target_class=a.target_class,diagnostics_dir=None if a.no_diagnostics else target); u.to_csv(o/'class_counts_timeseries_corrected.csv',index=False); write_long_uncertainty_csv(u,o/'class_counts_timeseries_corrected_long.csv'); visible=tuple(x.strip() for x in a.visible_classes.split(',')) if a.visible_classes else None; plot_timeseries(d,o/'class_counts_timeseries.html',a.stacked,a.log_y,visible_classes=visible,uncertainty_df=u)
if __name__=='__main__': main()
