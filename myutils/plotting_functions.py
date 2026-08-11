from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Local imports
from .general_functions import measureModelIR
from .model_functions import gov, customNomogram
from .nonlinearities import weib
from .spectral_templates import getNomogramParams

# Colors
aGS_color = '#ffa500' # orange. Was '#7300ff' purple
GS_color = '#32cd32' # lime green. Was #fe00fa' pink

# Fonts
plt.rcParams['font.family'] = 'Arial'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.size'] = 14

def sanityCheck():
    """Print a short message confirming this module loaded; used for notebook sanity checks."""
    print("plotting_functions.py Sanity Checked")


def _save_model_vs_data_pdf(fig, pdffile, crop):
    """Save ``fig`` to ``pdffile`` (PDF). Creates parent directories if needed."""
    if not pdffile:
        return
    path = Path(pdffile)
    path.parent.mkdir(parents=True, exist_ok=True)
    kwargs = {"format": "pdf", "bbox_inches": "tight"}
    if crop:
        kwargs["pad_inches"] = 0
    fig.savefig(path, **kwargs)


def _resolve_model_vs_data_tiff_path(pdffile, tifffile):
    """Return TIFF path, or None to skip. ``tifffile=False`` disables."""
    if tifffile is False:
        return None
    if tifffile not in (None, ""):
        return Path(tifffile)
    if pdffile:
        return Path("figures") / f"{Path(pdffile).stem}.tif"
    return None


def _save_model_vs_data_tiff(
    fig,
    tifffile,
    *,
    crop_mode,
    axs=None,
    data_time=None,
    model_time=None,
    plotdata=True,
    tiff_dpi=600,
):
    """Save raster TIFF under ``figures/`` with tight bbox; x-limits match data extent (normal mode)."""
    if not tifffile:
        return
    path = Path(tifffile)
    path.parent.mkdir(parents=True, exist_ok=True)
    save_kw = {
        "format": "tif",
        "dpi": tiff_dpi,
        "bbox_inches": "tight",
        "pad_inches": 0,
    }
    if crop_mode:
        fig.savefig(path, **save_kw)
        return

    ax0, ax1 = axs[0], axs[1]
    x0lim, x1lim = ax0.get_xlim(), ax1.get_xlim()
    if plotdata and len(data_time):
        t_lo = float(np.min(data_time))
        t_hi = float(np.max(data_time))
        if len(model_time):
            t_lo = min(t_lo, float(np.min(model_time)))
            t_hi = max(t_hi, float(np.max(model_time)))
    else:
        t_lo = float(np.min(model_time))
        t_hi = float(np.max(model_time))
    ax0.set_xlim(t_lo, t_hi)
    ax1.set_xlim(t_lo, t_hi)
    fig.savefig(path, **save_kw)
    ax0.set_xlim(x0lim)
    ax1.set_xlim(x1lim)


def _overlay_norm_factors(
    data_time,
    data_vals,
    model_time,
    prediction,
    model_norm_max_time_min=None,
):
    """
    Max-abs normalization factors for data vs model overlay.

    Uses the intersection of the data and model time supports (where both traces
    overlap). Factors are ``max(|values|)`` within that window only; callers
    apply them to the full arrays they plot.

    If ``model_norm_max_time_min`` is set, ``model_nf`` uses only model samples
    with ``model_time <= model_norm_max_time_min`` (same units as ``model_time``,
    i.e. minutes when ``model_time`` comes from ``xaxis / 60``).

    Returns
    -------
    data_nf : float
    model_nf : float
    """
    if model_time.size == 0:
        return 1.0, 1.0
    if data_time.size == 0:
        model_nf = float(np.nanmax(np.abs(prediction)))
        if not np.isfinite(model_nf) or model_nf <= 0:
            model_nf = 1.0
        return 1.0, model_nf

    dt_lo, dt_hi = float(np.min(data_time)), float(np.max(data_time))
    mt_lo, mt_hi = float(np.min(model_time)), float(np.max(model_time))
    t_lo = max(dt_lo, mt_lo)
    t_hi = min(dt_hi, mt_hi)
    if t_hi < t_lo:
        t_lo, t_hi = min(dt_lo, mt_lo), max(dt_hi, mt_hi)

    dw = (data_time >= t_lo) & (data_time <= t_hi)
    mw = (model_time >= t_lo) & (model_time <= t_hi)

    dv = data_vals[dw] if np.any(dw) else data_vals
    mw_model = mw
    if model_norm_max_time_min is not None:
        mw_cap = model_time <= float(model_norm_max_time_min)
        mw_model = mw & mw_cap
        if not np.any(mw_model):
            mw_model = mw
    pv = prediction[mw_model] if np.any(mw_model) else prediction

    data_nf = float(np.nanmax(np.abs(dv)))
    model_nf = float(np.nanmax(np.abs(pv)))
    if not np.isfinite(data_nf) or data_nf <= 0:
        data_nf = 1.0
    if not np.isfinite(model_nf) or model_nf <= 0:
        model_nf = 1.0
    return data_nf, model_nf


def plotModelVsData(
    prediction,
    data,
    xaxis,
    model,
    figsize=(10, 4),
    figfile="",
    pdffile="",
    tifffile=None,
    modelcolor=aGS_color,
    datacolor= "#4d4d4d",
    datatimebase=None,
    title="",
    scalebar=False,
    ylimit=None,
    xlimit=None,
    plotdata=True,
    modelalpha=0.5,
    dataalpha=0.8,
    crop=False,   # NEW
    model_line_width=2,
    tiff_dpi=600,
    model_norm_max_time_sec=None,
):
    """
    Plot melanopsin model prediction against experimental traces.

    In normal mode, builds a two-row figure: log10 light monitor (top) and normalized
    model current vs data (bottom). In ``crop`` mode, returns a single minimal axis
    suitable for figure panels (no axes decoration).

    Data and model are scaled by ``max(|value|)`` each, using only samples whose times
    lie in the intersection of the data time span and the model time span (after the
    model is cropped to ``model_time <= max(data_time)``). The same factors scale the
    full plotted traces. The model trace is then negated for voltage-clamp convention.

    Parameters
    ----------
    prediction : array_like
        Model output time series (e.g. current-related trace) aligned with ``xaxis``.
    data : array_like
        Experimental response; NaNs are dropped for plotting.
    xaxis : array_like
        Model time base in seconds (converted internally to minutes where noted).
    model : dict
        Must include ``lightIntensity`` (1D or use ``stimMonitor`` when 2D) for the light strip.
    figsize : tuple, optional
        Figure size in inches.
    figfile : str, optional
        If set, save figure to this path (PNG in normal mode, tight bbox in crop mode).
    pdffile : str or pathlib.Path, optional
        If set, save the same figure as a PDF (vector). Parent directories are created
        as needed. Saved after drawing; use with ``figfile`` if both formats are wanted.
    tifffile : str, pathlib.Path, or None, optional
        If set, save a TIFF to this path. If ``None`` and ``pdffile`` is set, saves
        ``figures/<pdf_stem>.tif`` with ``bbox_inches='tight'`` and ``pad_inches=0``;
        in normal mode the x-axis is set to the exact min/max of plotted data and model
        times for that export only. Pass ``False`` to disable TIFF when ``pdffile`` is set.
    tiff_dpi : int
        DPI for TIFF output (default 600).
    modelcolor, datacolor : str, optional
        Line colors for model and data.
    datatimebase : array_like or None
        Time base for data in seconds; if empty and ``plotdata``, inferred from ``xaxis`` span.
    title : str, optional
        Subplot title on the main axis.
    scalebar : bool
        If True, hide ticks and most spines on the main axis (scale-bar style).
    ylimit, xlimit : list or None
        Optional y limits for main axis; x limits for both rows.
    plotdata : bool
        Whether to plot experimental series.
    modelalpha, dataalpha : float
        Line alpha for model and data (defaults: model 0.8, data 0.5).
    crop : bool
        If True, single-axis cropped plot with normalized model and data only.
        Y limits use a small margin below (and a smaller margin above) the data
        range so rounded line caps are not clipped when saving with tight bbox.
    model_line_width : float
        Linewidth for the model trace.
    model_norm_max_time_sec : float or None
        If set, ``model_nf`` is ``max(|prediction|)`` using only times with
        ``xaxis <= model_norm_max_time_sec`` (seconds, same as stimulus ``timings``
        in ``setupModelRun``). Data scaling is unchanged.

    Returns
    -------
    matplotlib.figure.Figure
        The figure handle (``axs`` not returned; use ``fig.axes`` if needed).
    """
    # --------------------------------------------------
    # Input normalization
    # --------------------------------------------------
    data = np.asarray(data)
    prediction = np.asarray(prediction)

    if datatimebase is None:
        datatimebase = []
    if ylimit is None:
        ylimit = []
    if xlimit is None:
        xlimit = []

    _model_norm_max_time_min = (
        float(model_norm_max_time_sec) / 60.0
        if model_norm_max_time_sec is not None
        else None
    )

    # --------------------------------------------------
    # Construct data timebase
    # --------------------------------------------------
    if plotdata and len(datatimebase) == 0:
        datatimestep = np.max(xaxis) / len(data)
        datatimebase = np.arange(len(data)) * datatimestep

    not_nan_mask = ~np.isnan(data)
    data_time = datatimebase[not_nan_mask] / 60
    data_vals = data[not_nan_mask]

    # --------------------------------------------------
    # CROPPED MODE (single axis, no decoration)
    # --------------------------------------------------
    if crop:
        fig, ax = plt.subplots(figsize=figsize)

        # Model time crop (same x extent as before: model not past last data time)
        model_time = xaxis / 60
        model_mask = model_time <= np.max(data_time)
        model_time = model_time[model_mask]
        prediction = prediction[model_mask]

        data_nf, model_nf = _overlay_norm_factors(
            data_time,
            data_vals,
            model_time,
            prediction,
            model_norm_max_time_min=_model_norm_max_time_min,
        )

        data_vals_plot = np.array([], dtype=float)
        if plotdata:
            data_vals_plot = data_vals / data_nf
            ax.plot(
                data_time,
                data_vals_plot,
                color=datacolor,
                alpha=dataalpha,
            )

        # Negate model to match voltage-clamp sign convention
        norm_prediction = -prediction / model_nf

        ax.axhline(0, color="gray", dashes=[2, 2])

        ax.plot(
            model_time,
            norm_prediction,
            color=modelcolor,
            linewidth=model_line_width,
            alpha=modelalpha,
        )

        # Tight limits from actual plotted values
        x_all = np.concatenate([data_time, model_time])
        y_all = (
            np.concatenate([data_vals_plot, norm_prediction])
            if plotdata
            else norm_prediction
        )

        ax.set_xlim(np.min(x_all), np.max(x_all))
        y_min, y_max = float(np.min(y_all)), float(np.max(y_all))
        y_span = y_max - y_min
        if not np.isfinite(y_span) or y_span <= 0:
            y_span = 1.0
        # Tight bbox + exact ylim clips round line caps at extrema; small pad fixes it.
        pad_bottom = 0.025 * y_span
        pad_top = 0.01 * y_span
        ax.set_ylim(y_min - pad_bottom, y_max + pad_top)

        # Remove everything
        ax.set_axis_off()

        plt.tight_layout(pad=0)

        if figfile:
            plt.savefig(figfile, dpi=1000, bbox_inches="tight", pad_inches=0)
        _save_model_vs_data_pdf(fig, pdffile, crop=True)
        _tifpath = _resolve_model_vs_data_tiff_path(pdffile, tifffile)
        _save_model_vs_data_tiff(
            fig,
            _tifpath,
            crop_mode=True,
            plotdata=plotdata,
            tiff_dpi=tiff_dpi,
        )

        return fig

    # --------------------------------------------------
    # NORMAL MODE (original behavior)
    # --------------------------------------------------
    fig, axs = plt.subplots(
        2,
        1,
        figsize=figsize,
        gridspec_kw={"height_ratios": [0.1, 1], "hspace": 0.01},
    )

    # --- Light monitor ---
    if model["lightIntensity"].ndim == 1:
        log_light = np.log10(model["lightIntensity"] + 1)
    else:
        log_light = np.log10(model["stimMonitor"])

    light_time = xaxis / 60
    in_data_range = light_time <= np.max(data_time)

    axs[0].plot(
        light_time[in_data_range],
        log_light[in_data_range],
        color="k",
        linewidth=1.2,
    )

    ymin = np.nanmin(log_light[log_light > 0])
    ymax = np.nanmax(log_light)
    pad = 0.03 * (ymax - ymin)

    axs[0].set_ylim([ymin - pad, ymax + pad])

    axs[0].vlines(
        light_time[in_data_range],
        ymin=ymin,
        ymax=log_light[in_data_range],
        color="k",
        linewidth=0.5,
        linestyles="dotted",
        alpha=0.4,
    )

    for spine in ["top", "right", "left", "bottom"]:
        axs[0].spines[spine].set_visible(False)

    axs[0].set_xticks([])
    axs[0].set_yticks([])

    # --- Main plot ---
    axs[1].axhline(0, color="gray", dashes=[2, 2])

    model_time = xaxis / 60
    model_mask = model_time <= np.max(data_time)

    model_time = model_time[model_mask]
    prediction = prediction[model_mask]

    data_nf, model_nf = _overlay_norm_factors(
        data_time,
        data_vals,
        model_time,
        prediction,
        model_norm_max_time_min=_model_norm_max_time_min,
    )

    if plotdata:
        plot_data = np.asarray(data, dtype=float) / data_nf
        axs[1].plot(
            datatimebase / 60,
            plot_data,
            color=datacolor,
            alpha=dataalpha,
            label="Data",
        )

    norm_prediction = -prediction / model_nf

    axs[1].plot(
        model_time,
        norm_prediction,
        color=modelcolor,
        linewidth=model_line_width,
        alpha=modelalpha,
        label="Model",
    )

    if title:
        axs[1].set_title(title, y=1.06, pad=20)

    for spine in ["top", "right"]:
        axs[1].spines[spine].set_visible(False)

    if scalebar:
        for spine in ["top", "right", "bottom", "left"]:
            axs[1].spines[spine].set_visible(False)
        axs[1].set_xticks([])
        axs[1].set_yticks([])
    else:
        axs[1].set_xlabel("Time (min)")
        axs[1].set_ylabel("Current\n(norm.)")

    if ylimit:
        axs[1].set_ylim(ylimit)

    if xlimit:
        axs[0].set_xlim(xlimit)
        axs[1].set_xlim(xlimit)
    else:
        padding = 1
        axs[0].set_xlim([np.min(data_time) - padding, np.max(data_time) + padding])
        axs[1].set_xlim([np.min(data_time) - padding, np.max(data_time) + padding])

    axs[1].legend(bbox_to_anchor=(1.05, 1), loc="upper left")

    plt.tight_layout()
    # Leave room on the right for the outside legend (Tk canvas clips otherwise).
    fig.subplots_adjust(right=0.82)

    if figfile:
        plt.savefig(figfile, dpi=1000)
    _save_model_vs_data_pdf(fig, pdffile, crop=False)
    _tifpath = _resolve_model_vs_data_tiff_path(pdffile, tifffile)
    _save_model_vs_data_tiff(
        fig,
        _tifpath,
        crop_mode=False,
        axs=axs,
        data_time=data_time,
        model_time=model_time,
        plotdata=plotdata,
        tiff_dpi=tiff_dpi,
    )

    return fig

def plotIntensityResponseRelation(
    irdata,
    polarity,
    intensities,
    modelpredict,
    title,
    exampleidx,
    maxdata,
    maxmodel,
    figsize=(2.2, 2),
    figfile="",
    modelcolor=aGS_color,
    ax=None,
    ylabel="Response (Norm.)",
):
    """
    Plot one intensity–response panel: experimental traces vs model, log10 intensity on x-axis.

    Expects a dataframe-like ``irdata`` with columns ``Intensity`` and ``V*`` cellular measurements.
    One trace (``exampleidx``) is highlighted; others are faint. Model curve is normalized
    by ``maxmodel`` (scalar or per-row).

    Parameters
    ----------
    irdata : pandas.DataFrame
        Must include ``Intensity`` and columns whose names start with ``V``.
    polarity : float
        Sign applied to data (e.g. -1 for transient inward current).
    intensities : sequence
        Model intensities; values <= 0 are skipped on the x-axis.
    modelpredict : array_like
        Model responses at those intensities (same structure as ``maxmodel``).
    title : str
        Panel title.
    exampleidx : int
        Which ``V*`` column (1-based among V-columns) to emphasize.
    maxdata : dict-like
        Per-column maxima for normalizing experimental traces (keys match columns).
    maxmodel : float or sequence
        Normalization for ``modelpredict``.
    figsize : tuple, optional
        Used only if ``ax`` is None.
    figfile : str, optional
        If standalone, save figure to this path when provided.
    modelcolor : str, optional
        Model line color.
    ax : matplotlib.axes.Axes or None
        If provided, draw on this axis; otherwise create a new figure.
    ylabel : str, optional
        Y-axis label.

    Returns
    -------
    tuple
        If ``ax`` was provided: ``(alldata, dataintensity, norm_modelpredict, modelintensity)``.
        If a new figure was created: same four values plus ``fig`` as a fifth element.
    """
    # --------------------------------------------------
    # Figure / axis handling
    # --------------------------------------------------
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
        standalone_plot = True
    else:
        standalone_plot = False

    # --------------------------------------------------
    # Prepare data
    # --------------------------------------------------
    dataintensity = np.log10(irdata["Intensity"])

    alldata = []
    dataintensity_example = None
    dataresponse_example = None

    # --------------------------------------------------
    # Plot experimental data
    # --------------------------------------------------

    trace_idx = 1

    for column in irdata.columns:
        if not column.startswith("V"):
            continue

        dataresponse = polarity * irdata[column] / maxdata[column]
        alldata.append(dataresponse)

        if trace_idx == exampleidx:
            dataintensity_example = dataintensity
            dataresponse_example = dataresponse
        else:
            ax.plot(
                dataintensity,
                dataresponse,
                color="#6e6e6e",
                linewidth=0.72,
                zorder=1
            )

        trace_idx += 1

    # Highlight example trace if requested
    if dataintensity_example is not None:
        ax.plot(
            dataintensity_example,
            dataresponse_example,
            color="#4d4d4d",
            marker="o",
            markersize=5,
            markerfacecolor="w",
            linewidth=2,
            zorder=1
        )

    # --------------------------------------------------
    # Plot model predictions
    # --------------------------------------------------
    modelintensity = np.log10([i for i in intensities if i > 0])

    if np.ndim(maxmodel) == 0:
        norm_modelpredict = modelpredict / maxmodel
        ax.plot(
            modelintensity,
            norm_modelpredict,
            color=modelcolor,
            linewidth=2,
            marker=".",
            label="Model",
            zorder=3
        )
    else:
        norm_modelpredict = []
        for i, mm in enumerate(maxmodel):
            nm = modelpredict[i] / mm
            norm_modelpredict.append(nm)
            ax.plot(
                modelintensity,
                nm,
                color=modelcolor[i],
                linewidth=2,
                label="Model",
                zorder=3
            )

    # --------------------------------------------------
    # Axis formatting
    # --------------------------------------------------
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Y limits
    ymax = np.max(norm_modelpredict)
    ax.set_ylim([-0.01, 1.01 if ymax > 0.2 else 0.2])

    # X ticks: integers strictly within data range
    xmin, xmax = ax.get_xlim()

    int_ticks = np.arange(
        int(np.ceil(xmin)),
        int(np.floor(xmax)) + 1,
        1,
    )

    ax.set_xticks(int_ticks)

    ax.set_xlabel("Log photons")
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    # --------------------------------------------------
    # Finalize
    # --------------------------------------------------
    if standalone_plot:
        if figfile:
            plt.savefig(figfile)
        plt.show()
        return (
            alldata,
            dataintensity,
            norm_modelpredict,
            modelintensity,
            fig,
        )

    return alldata, dataintensity, norm_modelpredict, modelintensity

def plotIntensityResponseRelations(
    irt,
    irp,
    intensities,
    transient,
    persistent,
    replength,
    maxdata,
    exampleidx=-1,
    modelcolor=aGS_color,
    transientb=None,
    persistentb=None,
    modelbcolor=GS_color,
    grid2x2=False,
):
    """
    Plot four intensity–response panels: transient and persistent, split across two stimulus sets.

    Uses :func:`plotIntensityResponseRelation` for each panel. Optionally overlays a second
    model condition (e.g. background) via ``transientb`` / ``persistentb``.

    Parameters
    ----------
    irt, irp : pandas.DataFrame
        Intensity–response data for transient (negative polarity) and persistent (positive)
        responses; filtered by ``step`` vs ``replength``.
    intensities : sequence
        Full intensity list; split at ``replength`` for the two sets.
    transient, persistent : array_like
        Model predictions (vector or 2D with one row per condition).
    replength : int
        Number of stimuli in the first set; remainder is the second set.
    maxdata : dict-like
        Normalization for experimental traces (passed through to each panel).
    exampleidx : int, optional
        Which example trace to highlight in each panel (see :func:`plotIntensityResponseRelation`).
    modelcolor, modelbcolor : str, optional
        Colors for primary and secondary model overlays.
    transientb, persistentb : array_like or list, optional
        Optional secondary model curves for the transient and persistent panels.
    grid2x2 : bool
        If True, use a 2×2 layout with shared x-axis and square aspects; if False, 1×4.

    Returns
    -------
    tuple
        ``(alldata1t, alldata2t, alldata1p, alldata2p, fig)`` — experimental list aggregates
        per panel and the figure handle.
    """
    if transientb is None:
        transientb = []
    if persistentb is None:
        persistentb = []

    # --------------------------------------------------
    # Figure / axes creation
    # --------------------------------------------------
    if grid2x2:
        fig, axs = plt.subplots(
            2, 2,
            figsize=(6, 6),
            sharex=True
        )
        axs = axs.flatten()
    else:
        fig, axs = plt.subplots(1, 4, figsize=(10, 3))
        axs = axs.flatten()

    # --------------------------------------------------
    # Data shape handling
    # --------------------------------------------------
    is_vector = len(np.shape(transient)) == 1
    maxmodel = np.max(transient) if is_vector else np.max(transient, axis=1)

    # Normalize persistentb to global max of transientb
    if transientb is not None and len(transientb) > 0:
        if is_vector:
            max_transientb = np.max(transientb)
        else:
            max_transientb = np.max(transientb, axis=1)
    else:
        max_transientb = 1

    # --------------------------------------------------
    # Helpers
    # --------------------------------------------------
    def overlay_secondary(ax, secondary, intensities_subset, normalize_to):
        if secondary is None or len(secondary) == 0:
            return

        modelintensity = np.log10([i for i in intensities_subset if i > 0])

        if is_vector:
            ax.plot(
                modelintensity,
                secondary / normalize_to,
                color=modelbcolor,
                linewidth=2,
                alpha=1,
                marker="."
            )
        else:
            for i, row in enumerate(secondary):
                ax.plot(
                    modelintensity,
                    row / normalize_to[i],
                    color=modelbcolor,
                    linewidth=2,
                    alpha=1,
                    marker=".",
                    zorder=2
                )

    def plot_panel(
      irdata,
      modelpredict,
      intensities_subset,
      polarity,
      ax,
      secondary=None,
      title=""
  ):
      # Plot secondary (background) first
      overlay_secondary(ax, secondary, intensities_subset, max_transientb)

      # Plot primary on top
      alldata, _, _, _ = plotIntensityResponseRelation(
          irdata,
          polarity,
          intensities_subset,
          modelpredict,
          title,
          exampleidx,
          maxdata,
          maxmodel,
          modelcolor=modelcolor,
          ax=ax,
          ylabel="Peak Response"
      )

      if grid2x2:
          ax.set_box_aspect(1)
          ax.set_xlabel("")
          ax.set_ylabel("")

      return alldata


    # --------------------------------------------------
    # Data splits
    # --------------------------------------------------
    intensities1 = intensities[:replength]
    intensities2 = intensities[replength:]

    transient1 = transient[:replength] if is_vector else transient[:, :replength]
    transient2 = transient[replength:] if is_vector else transient[:, replength:]

    persistent1 = persistent[:replength] if is_vector else persistent[:, :replength]
    persistent2 = persistent[replength:] if is_vector else persistent[:, replength:]

    transientb1 = transientb[:replength] if is_vector else transientb[:, :replength]
    transientb2 = transientb[replength:] if is_vector else transientb[:, replength:]

    persistentb1 = persistentb[:replength] if is_vector else persistentb[:, :replength]
    persistentb2 = persistentb[replength:] if is_vector else persistentb[:, replength:]

    # --------------------------------------------------
    # Panels
    # --------------------------------------------------
    alldata1t = plot_panel(
        irt[irt.step < replength + 1],
        transient1,
        intensities1,
        -1,
        axs[0],
        secondary=transientb1,
        title="Set 1"
    )

    alldata2t = plot_panel(
        irt[irt.step > replength],
        transient2,
        intensities2,
        -1,
        axs[1],
        secondary=transientb2,
        title="Set 2"
    )

    alldata1p = plot_panel(
        irp[irp.step < replength + 1],
        persistent1,
        intensities1,
        1,
        axs[2],
        secondary=persistentb1,
        title="" if grid2x2 else "Set 1"
    )

    alldata2p = plot_panel(
        irp[irp.step > replength],
        persistent2,
        intensities2,
        1,
        axs[3],
        secondary=persistentb2,
        title="" if grid2x2 else "Set 2"
    )

    # --------------------------------------------------
    # Axis cleanup for 2x2 layout
    # --------------------------------------------------
    if grid2x2:
        # Remove x tick labels on top row
        for ax in axs[:2]:
            ax.tick_params(labelbottom=False)

        # Keep x tick labels on bottom row
        for ax in axs[2:]:
            ax.tick_params(labelbottom=True)

        # Remove y tick labels on right column
        for ax in (axs[1], axs[3]):
            ax.tick_params(labelleft=False)

        # Shared labels
        fig.supxlabel("Log Intensity", y=0.04)
        axs[0].set_ylabel("Peak Response")
        axs[2].set_ylabel("Peak Response")

        # Thicker axis spines and tick marks


    plt.tight_layout()
    plt.show()

    return alldata1t, alldata2t, alldata1p, alldata2p, fig

def plotConditioning(model, stim, figsize=(4, 2), xlimit=None):
    """
    Plot 480:440 spectral sensitivity ratio vs time of 560 nm conditioning.

    Extracts transient responses from ``model['currentGlobalGain']`` at stimulus onsets
    in ``stim``, compares 440 nm and 480 nm probe responses, and plots the ratio timecourse
    with reference lines for pure R vs pure E templates.

    Parameters
    ----------
    model : dict
        Output of the melanopsin simulation; uses ``currentGlobalGain``, ``rate``, and
        timing-compatible fields.
    stim : dict
        Stimulus description with ``timings``, ``interstim``, ``wlen`` (wavelength per pulse).
    figsize : tuple, optional
        Figure size in inches (default ``(4, 2)``).
    xlimit : list or None, optional
        X-axis limits ``[xmin, xmax]``. Auto-scaled if *None*.

    Returns
    -------
    matplotlib.figure.Figure
        Current pyplot figure after plotting.
    """

    fig = plt.figure(figsize=figsize)

    start_times = np.array(stim['timings'][::2])
    search_window = 10
    interstim = stim['interstim']
    transient, persistent, prestimbaselines = measureModelIR(
        model['currentGlobalGain'],
        start_times,
        search_window,
        interstim,
        model['rate'],
    )

    start_times = np.array(stim['timings'][::2]) / 60  # convert sec to min
    start_time_560 = start_times[np.array(stim['wlen']) == 560]
    start_times = start_times - start_time_560
    stims_440 = np.array(stim['wlen']) == 440
    stims_480 = np.array(stim['wlen']) == 480

    responses_440 = np.array(transient)[stims_440]
    responses_480 = np.array(transient)[stims_480]

    # Theoretical 480:440 sensitivity ratios for pure R and pure E
    params = model['params']
    wl_probe = np.array([440.0, 480.0])
    nomoR = customNomogram(wl_probe, getNomogramParams(params, 'R'))
    nomoE = customNomogram(wl_probe, getNomogramParams(params, 'E'))
    Ract_480_440 = float(nomoR[1] / nomoR[0])
    Eact_480_440 = float(nomoE[1] / nomoE[0])

    ratio = [r480 / r440 for r440, r480 in zip(responses_440, responses_480)]
    ratio_times = start_times[stims_440]
    ratio_times_filt = np.delete(ratio_times, 10)
    ratio_filt = np.delete(ratio, 10)

    plt.plot(ratio_times_filt, ratio_filt, '.k', label='Model')
    plt.axvline(0, ls='--', c='gray')

    plt.errorbar(-10, 1.15, 0.05, fmt='dk', label='Emanuel 2015 (dark)')
    plt.errorbar(5, 0.78, 0.05, fmt='dr', label='Emanuel 2015 (conditioned)')

    plt.ylabel('Sensitivity Ratio\n(480:440)')
    plt.xlabel('Time (min.)')

    mean_before = np.mean(ratio[:10])
    mean_after = np.mean(ratio[11:15])

    plt.axhline(mean_before, ls='--', c='gray')
    plt.axhline(Ract_480_440, ls='dotted', c='k', label='R')
    plt.axhline(Eact_480_440, ls='dotted', c='r', label='E')

    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    max_time = np.floor(np.max(ratio_times_filt) / 10) * 10
    plt.ylim([0.6, 1.25])
    if xlimit is None:
        plt.xlim([-15, max_time])
    else:
        plt.xlim(xlimit)
    plt.xticks(np.arange(-10, max_time, 10))

    print('480:440 ratio before conditioning: ' + str(mean_before))
    print('480:440 ratio after conditioning: ' + str(mean_after))
    print('Alan measured 1.15 before and 0.78 after')

    return fig

def overlayModelActivity(models, labels=None, alpha=0.7):
    """
    Overlay ``currentGlobalGain`` from several model runs on one time base.

    Parameters
    ----------
    models : sequence of dict
        Each dict must include ``xaxis`` and ``currentGlobalGain`` (as from the simulator).
    labels : sequence of str or None
        Legend labels; if None, uses empty strings (one per model).
    alpha : float
        Line transparency.

    Returns
    -------
    matplotlib.figure.Figure
        The figure handle from ``plt.figure()``.
    """
    h = plt.figure()

    if labels == None:
        labels = "" * len(models)

    for model,label in zip(models, labels):
        timebase = model['xaxis']
        activity = model['currentGlobalGain']
        plt.plot(timebase, activity, label=label, alpha=alpha)

    plt.legend()

    plt.xlabel('Time (s)')
    plt.ylabel('Activity')

    return h 

def plotSpectralTemplates(param_list, labels=None, wlen=None, colors=None, logy=False):
    """
    Plot spectral templates generated from customNomogram() or gov().

    Parameters
    ----------
    param_list : list of dict
        List of parameter dictionaries.
    labels : list of str, optional
        Labels for each template.
    wlen : ndarray, optional
        Wavelength vector. If None, defaults to 350–650 nm.
    colors : list, optional
        Colors for each template.
    logy : bool, optional
        If True, plot sensitivity on a log y-axis.

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes
    """

    if wlen is None:
        wlen = np.linspace(350, 650, 1000)

    if labels is None:
        labels = [f"Template {i+1}" for i in range(len(param_list))]

    if colors is None:
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color'][:len(param_list)]

    fig, ax = plt.subplots()

    for params, label, color in zip(param_list, labels, colors):

        spec_type = params.get('spectral_sensitivity', 'Custom')

        if spec_type == 'Govardovskii':
            lambdamax = params['true_lambda_max']
            nomo = gov(wlen, lambdamax)

        else:
            nomo = customNomogram(wlen, params)

        if logy:
            ax.semilogy(wlen, nomo, label=label, color=color)
        else:
            ax.plot(wlen, nomo, label=label, color=color)

    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Normalized spectral sensitivity")
    ax.set_title("Spectral Templates")
    ax.legend()
    ax.set_ylim(0, 1.05)

    return fig, ax


def sigmoidVectorized(x, k, x0):
    """
    Description
    -----------
    3pl sigmoid for adding nonlinearity to the model. Full description found at https://en.wikipedia.org/wiki/Logistic_function

    Parameters
    ----------
    x : array_like
        Input values (scalar or array).
    k : float
        steepness
    x0 : float
        "Midpoint" of the curve, or its inflection point, or the value halfway to maximum.

    Returns
    -------
    Sigmoidal value at point x (float).
    """
    a = k * (x - x0)
    return np.where(
        a < 0,
        np.exp(a) / (1 + np.exp(a)),         # if a < 0
        1 / (1 + np.exp(-a))                     # if a >= 0
    )

def plotSigmoid(k = 40, x0 = 8e-4):
    """
    Plot the activity-feedback nonlinearity: ``1 - sigmoidVectorized(x, k, x0)`` vs ``x``.

    Parameters
    ----------
    k : float
        Steepness of the logistic.
    x0 : float
        Midpoint / inflection parameter.

    Returns
    -------
    None
        Displays the figure via ``plt.show()``.
    """
    x = np.linspace(-0.5, 1, 10000)
    y = 1-sigmoidVectorized(x, k, x0)
    plt.figure(figsize=(8, 5))
    plt.plot(x, y)
    plt.title(f'Plot with k={k}, x0={x0}')
    plt.ylim([0,1.1])
    plt.xlabel('Feedback from Activity')
    plt.ylabel('Global Gain')
    plt.grid(True)
    plt.show()


def plotWeib(k=1.215, x0=0.381):
    """
    Plot the M′-driven feedback shape: Weibull CDF from :func:`~myutils.nonlinearities.weib`.

    Parameters
    ----------
    k, x0 : float
        Weibull shape and scale passed to ``weib``.

    Returns
    -------
    None
        Uses pyplot state; does not call ``plt.show()`` (display depends on the environment).
    """
    x = np.linspace(0, 3, 100000)
    y = weib(x, k, x0)
    plt.figure(figsize=(8,5))
    plt.plot(x, y)
    plt.title(f'Plot with, k={k}, x0={x0}')
    plt.ylim([0, 1])
    plt.xlabel('Feedback from M prime')
    plt.ylabel('M gain')
    plt.grid(True)


def plotModelRun(model, figsize=[10,10], title = '', figfile='',
                 plotcurrent = False, Mcolor='#663695', MprimeColor='#D28ABB', Rcolor='k', RprimeColor='#808080', Ecolor='#EC1E24',
                 EprimeColor='#F37F72', Ocolor='#6a2a00', globalgaincolor='#2ca02c', Mgaincolor='#f05928', activity_ylimit=None,
                 pigment_fraction_ylim_fixed=True, show=True):
    """
    Multi-panel summary of a melanopsin simulation: light, activity, pigment states, optional gains.

    Time is shown in minutes. Panels typically include log10 light, ``currentGlobalGain``,
    M/M′, R/R′, E/E′, O, and (if present) normalized global and M gain scalars.

    Parameters
    ----------
    model : dict
        Output dict from the simulator (e.g. :func:`~myutils.model_functions.predictMelanopsin`):
        ``xaxis``, ``currentGlobalGain``, ``currentOut``, pigment fractions ``M``, ``Mprime``,
        ``R``, ``Rprime``, ``E``, ``Eprime``, ``O``, ``lightIntensity`` or ``stimMonitor``,
        and optionally ``gg_scalar``, ``mg_scalar``.
    figsize : list of float, optional
        ``[width, height]`` in inches.
    title : str, optional
        Figure suptitle.
    figfile : str, optional
        If non-empty, save as PDF to this path.
    show : bool, optional
        If True (default), call ``plt.show()`` when ``figfile`` is empty. If False, skip
        interactive display so the figure can be embedded (e.g. TkAgg canvas).
    plotcurrent : bool
        If True, overlay ``currentOut`` on the activity panel.
    Mcolor, MprimeColor, Rcolor, RprimeColor, Ecolor, EprimeColor, Ocolor : str
        Colors for pigment state traces.
    globalgaincolor, Mgaincolor : str
        Colors for gain traces in the bottom panel.
    activity_ylimit : tuple or None
        Fixed y-limits for the activity panel; if None, auto-scales to max activity.
    pigment_fraction_ylim_fixed : bool
        If True, pigment panels use y in ``[-0.1, 1.1]`` with y-axis ticks at 0, 0.5, and 1.

    Notes
    -----
    Default line widths match Figure 5d panels: light 1.0, activity and pigment traces 1.4,
    optional ``currentOut`` overlay 1.0 dashed, gain traces 1.4 when present.

    Returns
    -------
    fig : matplotlib.figure.Figure
    axs : ndarray of matplotlib.axes.Axes
        Vertical stack of subplots (length 6 or 7 depending on gain data).
    """
    # Unpack results from the model run
    xaxis = model['xaxis'] / 60 # Convert time to minutes
    currentGlobalGain = model['currentGlobalGain']
    current = model['currentOut']
    M = model['M']
    Mprime = model['Mprime']
    R = model['R']
    Rprime = model['Rprime']
    E = model['E']
    Eprime = model['Eprime']
    O = model['O']
    gg = model['gg']
    lightIntensity = model['lightIntensity']
    mg = model['mg']
    gg = model['gg']
    mg_scalar = model['mg_scalar']
    gg_scalar = model['gg_scalar']

    # Ensure xaxis is the correct length
    if len(xaxis) < len(currentGlobalGain):
        xaxis = np.append(xaxis, xaxis[-1] + np.median(np.diff(xaxis)))
        print('Added an element to xaxis to match length of melanopsin activity')
    elif len(xaxis) > len(currentGlobalGain):
        xaxis = xaxis[:-1]
        print('Removed an element from xaxis to match length of melanopsin activity')

    # Check if gain control was turned on during the model run
    if len(gg) > 0 and len(mg) > 0:
        nsubplots = 7
    else:
        nsubplots = 6

    fig, axs = plt.subplots(nsubplots, 1, figsize=(figsize[0], figsize[1]), gridspec_kw={'height_ratios':[0.5] + (nsubplots - 1) * [1]})

    if len(title) > 0:
        fig.suptitle(title, fontsize=16)
    else:
        pass

    if lightIntensity.ndim == 1:
        logLightIntensity = np.log10(np.maximum(lightIntensity, 0) + 1)
    else:
        stim_monitor = np.asarray(model['stimMonitor'])
        logLightIntensity = np.log10(np.maximum(stim_monitor, 0) + 1)

    _lw_light = 1.0
    _lw_activity = 1.4
    _lw_activity_current = 1.0
    _lw_pigment = 1.4
    _lw_O = 1.4
    _lw_gain = 1.4
    _lw_axhline = 1.0

    axs[0].plot(xaxis, logLightIntensity, color='k', linewidth=_lw_light)
    axs[0].set_ylabel('Light')
    # Scale y-axis to data range so intensity variation is visible (no fixed baseline at 0)
    finite_light = logLightIntensity[np.isfinite(logLightIntensity)]
    if finite_light.size == 0:
        ylo, yhi = 0.0, 1.0
    else:
        ylo = np.min(finite_light)
        yhi = np.max(finite_light)
    span = yhi - ylo
    if span < 1e-10:
        pad = max(1.0, abs(ylo) * 0.1)
        ylo, yhi = ylo - pad, yhi + pad
    else:
        pad = 0.1 * span
        ylo, yhi = ylo - pad, yhi + pad
    axs[0].set_ylim([ylo, yhi])

    axs[1].axhline(y=0, color='gray', linewidth=_lw_axhline, dashes=[5, 5])
    axs[1].plot(xaxis, currentGlobalGain, color='k', linewidth=_lw_activity)
    axs[1].set_ylabel('Activity')
    if plotcurrent:
        axs[1].plot(xaxis, current, label='Current (GS)', color='grey',
                    linewidth=_lw_activity_current, linestyle='--')
    if activity_ylimit is None:
        try:
            axs[1].set_ylim([0, 1.2*np.max(currentGlobalGain)])
        except:
            print('Failed to set axis limits')
    else:
        axs[1].set_ylim(activity_ylimit)

    axs[2].plot(xaxis, M, label='M', color=Mcolor, linewidth=_lw_pigment)
    axs[2].plot(xaxis, Mprime, label='M\'', color=MprimeColor, linewidth=_lw_pigment)

    axs[3].plot(xaxis, R, label='R', color=Rcolor, linewidth=_lw_pigment)
    axs[3].plot(xaxis, Rprime, label='R\'', color=RprimeColor, linewidth=_lw_pigment)

    axs[4].plot(xaxis, E, label='E', color=Ecolor, linewidth=_lw_pigment)
    axs[4].plot(xaxis, Eprime, label='E\'', color=EprimeColor, linewidth=_lw_pigment)

    axs[4].set_ylabel('Pigment\nFraction')

    axs[5].plot(xaxis, O, label='O', color=Ocolor, linewidth=_lw_O)

    if len(gg_scalar) > 0 and len(mg_scalar) > 0:
        axs[6].plot(xaxis, gg_scalar / np.max(gg_scalar), label='Global Gain', color=globalgaincolor,
                    linewidth=_lw_gain)  # normalized gains
        axs[6].plot(xaxis, mg_scalar / np.max(mg_scalar), label='M Gain', color=Mgaincolor,
                    linewidth=_lw_gain)
        axs[6].set_ylabel('Gain')
        axs[6].axhline(y=0.00, color='gray', linestyle='--', linewidth=_lw_axhline, alpha=1)

    for a in range(len(axs)):
        if a > 1:
            axs[a].legend(loc='upper left', bbox_to_anchor=(1, 1))

        if a < nsubplots - 1:
            axs[a].set_xticks([])
        axs[a].spines['top'].set_visible(False)
        axs[a].spines['right'].set_visible(False)
        if a > 0:
            axs[a].spines['bottom'].set_visible(False)
        axs[a].set_xlim([0,np.max(xaxis)])
        # Pigment fraction panels (M/M', R/R', E/E', O) are axes 2, 3, 4, 5
        if 2 <= a <= 5:
            if pigment_fraction_ylim_fixed:
                axs[a].set_ylim([-0.1, 1.1])
                axs[a].set_yticks([0.0, 0.5, 1.0])
                axs[a].set_yticklabels(["0", "0.5", "1"])

    axs[-1].spines['bottom'].set_visible(True)
    axs[0].spines['bottom'].set_visible(True)

    axs[-1].set_xlabel('Time (min.)')
    # Leave room on the right for legends (bbox_to_anchor=(1, 1) on pigment / gain panels).
    fig.subplots_adjust(
        left=0.09, right=0.72, top=0.93, bottom=0.07, hspace=0.14
    )

    if len(figfile) > 0:
        plt.savefig(figfile, format = 'pdf')
    elif show:
        plt.show()

    return fig, axs