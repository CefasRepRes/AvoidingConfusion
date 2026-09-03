#!/usr/bin/env python3
from itertools import cycle
import numpy as np, plotly.graph_objects as go, plotly.offline as pyo
from plotly.colors import qualitative,hex_to_rgb
META={'timestamp','run_name','blob_name'}
def plot_timeseries(df,path,stacked,log_y,y_min=None,y_max=None,visible_classes=None,uncertainty_df=None):
 fig=go.Figure(); colours=cycle(qualitative.Plotly); defaults=set(visible_classes or [c for c in df if c not in META])
 for c in df[[x for x in df if x not in META]].sum().sort_values(ascending=False).index:
  col=next(colours); vis=True if c in defaults else 'legendonly'; lo=f'{c}_corrected_lower'; hi=f'{c}_corrected_upper'; med=f'{c}_corrected_median'
  if uncertainty_df is not None and lo in uncertainty_df and hi in uncertainty_df:
   rgb=hex_to_rgb(col); fig.add_trace(go.Scatter(x=df.timestamp,y=uncertainty_df[hi],mode='lines',line={'width':0},visible=vis,showlegend=False,legendgroup=c)); fig.add_trace(go.Scatter(x=df.timestamp,y=uncertainty_df[lo],mode='lines',line={'width':0},fill='tonexty',fillcolor=f'rgba({rgb[0]},{rgb[1]},{rgb[2]},.2)',visible=vis,showlegend=False,legendgroup=c))
  fig.add_trace(go.Scatter(x=df.timestamp,y=df[c],mode='lines',name=f'{c} observed prediction',line={'color':col,'dash':'dot'},visible=vis,legendgroup=c))
  fig.add_trace(go.Scatter(x=df.timestamp,y=uncertainty_df[med] if uncertainty_df is not None else df[c],mode='lines+markers',name=c,line={'color':col},visible=vis,legendgroup=c,stackgroup='one' if stacked else None))
 fig.update_layout(title='Class counts through time across all runs',xaxis_title='Timestamp',yaxis_title='Count, log scale' if log_y else 'Count',hovermode='x unified',height=850); fig.update_xaxes(rangeslider_visible=True); fig.update_yaxes(type='log' if log_y else 'linear')
 pyo.plot(fig,filename=str(path),auto_open=False,include_plotlyjs='cdn')
