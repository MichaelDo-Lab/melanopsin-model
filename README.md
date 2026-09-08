# Melanopsin Model

Light regulates the circadian clock, mood, sleep, development, and vision. Underlying all of these functions are neural signals originating from melanopsin cells in the eye. Here, we present a model that predicts electrical signals called photocurrents that melanopsin cells generate in response to light. 

For an in-depth explanation of the biological context and the mechanisms instantiated in the model, see the open-access preprint: Nguyen and Caval-Holme et al. 2026.

The model runs through an app with a graphical user interface (GUI), an interactive tutorial notebook, or from source code. The GUI allows users to interact with the model immediately, without having to think about underlying implementation details. The tutorial notebook provides a more technical overview of stimulus creation, model parameter tuning, and ways to interact directly with the source code.

## Releases

This repository is released in two tracks:

- **bioRxiv preprint:** tag [`v0.1.0-bioRxiv2026`](https://github.com/MichaelDo-Lab/melanopsin-model/releases/tag/v0.1.0-bioRxiv2026) freezes the code and Windows executable that match Nguyen and Caval-Holme et al. 2026 (bioRxiv). Use this snapshot to reproduce the preprint.
- **Continuously updated:** the default branch (`main` or the current development branch) and the **Latest** release on the [Releases](https://github.com/MichaelDo-Lab/melanopsin-model/releases) page receive ongoing fixes and improvements. 

## Run the app (no Python required)

The GUI ships as a single Windows executable.

1. Download this repository or clone it using Git:

   ```bash
   git clone https://github.com/MichaelDo-Lab/melanopsin-model.git
   ```

2. Download `MelanopsinModel-v<version>.exe` from the repository's
   **Releases** page and place it in the downloaded/cloned folder.

3. Double-click the `.exe`.

The executable reads the spectral assets in `data/` and writes results to
`outputs/`. **Keep the `.exe` inside the cloned repository** so it can find
`data/` and `myutils/`. (It searches upward from its own location for those
folders, so placing it in the repository root or in `dist/` both work.)

No Python installation is required to run the executable.

## Run the tutorial notebook

The interactive tutorial is [`Melanopsin_Model_Tutorial.ipynb`](Melanopsin_Model_Tutorial.ipynb).

### Google Colab

Recommended if you want to try the model without installing Python.

1. Open [Google Colab](https://colab.research.google.com/).
2. Upload the notebook (`File > Upload notebook`), or open it from GitHub (`File > Open notebook > GitHub`) if the repository is available to Colab.
3. Run the first code cell under **Initialization**, or use **Runtime > Run all**. That cell detects Colab, clones the repository when needed, installs dependencies from `requirements.txt`, and loads the model.

Colab sessions are temporary—download any plots or data you want to keep.

### Local Jupyter

Requires Python 3.10.x. From the repository root:

```bash
git clone https://github.com/MichaelDo-Lab/melanopsin-model.git
cd melanopsin-model
pip install -r requirements.txt
pip install jupyter
jupyter notebook Melanopsin_Model_Tutorial.ipynb
```

You can also open the same `.ipynb` in a programming environment (VS Code/Cursor, Pycharm, JupyterLab, Spyder...) with a Python 3.10 kernel. Run the **Initialization** cell with the notebook in the repository root (the notebook will try to move there if needed).

## Run from source (for development)

Requires Python 3.10.x.

```bash
pip install -r requirements.txt
python run_melanopsin_gui.py
```

You can also launch the GUI as a module:

```bash
python -m myutils.melanopsin_gui
```

## Repository layout

Files and folders you will use:

| Path | Contents |
| --- | --- |
| `Melanopsin_Model_Tutorial.ipynb` | The interactive tutorial notebook (Google Colab or local Jupyter). |
| `data/` | Container for manuscript assets and the user library. |
| `data/manuscript/` | Container for manuscript recordings, spectra, and irradiance-response measurements. |
| `data/manuscript/recordings/` | Exemplar electrophysiological recording. |
| `data/manuscript/spectra/` | Measured light spectra (440 nm, 560 nm, Xenon, ocular filtering, and related assets). |
| `data/manuscript/irradiance response measurements/` | Recorded transient and persistent responses across light intensities, used to compare predictions against experiments. |
| `data/user_library/` | Spectra, light stimuli, and parameter presets you save in the app. |
| `outputs/` | Where your results are saved: figures in `images/` and exported prediction data in `predictions/`. This folder appears the first time you save something. |
| `run_melanopsin_gui.py` | Starts the GUI when you run `python run_melanopsin_gui.py`. |
| `requirements.txt` | The list of Python packages the model needs. |

Files for developers:

| Path | Contents |
| --- | --- |
| `myutils/` | The model code, plotting functions, and the desktop GUI (`melanopsin_gui.py`). |
| `myutils/app_paths.py` | Works out where `data/` and `outputs/` live, both when running from source and from the packaged app. |
| `myutils/_version.py` | The version number, which sets the executable's file name and window title. |
| `packaging/` | Scripts and settings for building the Windows executable (see [`packaging/README.md`](packaging/README.md)). |
| `.github/workflows/` | Automation that builds the executable and attaches it to a GitHub Release when a version tag such as `v0.2.0` is pushed. |
| `LICENSE` | Terms of use. |

## Building the executable

See [`packaging/README.md`](packaging/README.md). In short, from the repository
root on Windows:

```powershell
powershell -ExecutionPolicy Bypass -File packaging\build_windows.ps1 -Clean
```

This produces `dist/MelanopsinModel-v<version>.exe`.

# The GUI

A graphical user interface (GUI) allows users to create or import light spectra,
build light stimulus protocols, and export model predictions. 

## The Main Page: Selecting a Light Stimulus and Running the Model 

From the main window, the user can select light stimulus protocols and run the model. 

Upon first loading the GUI, the only available stimuli are those that were described in the manuscript.

At the top of the main window is a tool bar with elements that allows
the user to import/export data, configure model runs, create custom
spectra, and build light stimulus protocols. The exact functionalities
of this toolbar will be expanded upon in later sections.

In order to run a model prediction, follow these steps:

1. Click the stimulus bar at the top of the main window
   (upon launch it reads "Select stimulus to run..."). A drop-down
   list of all stimuli in the stimulus library opens in place.

2. Select a stimulus from the drop-down. The chosen stimulus appears
   in the **Stimulus:** bar.

3. Click the green **Run model** button to run the model.

While the model runs, a progress bar and a status message to the right of the buttons
report how far along the run is. The red **Stop** button cancels the current run. If an
earlier run finished, its results stay available for export.

## File

This menu allows the user to export data from a model run and import other
data such as stimuli and spectra.

### Export Model Run

The functionality of **Export Model Run** is only available after fully running a model prediction
and having the model plots shown on the **Main Page**. Once a model has successfully run
the user can export the prediction data by

1. Selecting **File** then **Export Model Run** (or **Ctrl + S** on the **Main Page**).

2. Activating **Export Model Run** will open a pop-up menu where the user 
   can select where the model data will be saved. By default it opens
   in `outputs/predictions`.

Data from a model prediction can be saved as a `.csv` table or as a `.npz`
compressed NumPy archive. Each column of the `.csv` is a model output
and each row is the value of that output at one time step.
Export includes only the time-series outputs (not the model's parameters in the
`params` dictionary). To keep the parameter values used for a run, save
them separately with **Save Parameters** in **Configure** →
**Model parameters…** (see below).

### Import Spectrum

**Import Spectrum** allows the user to import a new spectrum into the GUI.
These spectral files have specific requirements (see below).

1. Select **File** then **Import Spectrum**.

2. Choose a data file containing wavelength and intensity columns
   (``.csv``, ``.xlsx``, ``.xls``, ``.tsv``, or ``.txt``). Named columns
   such as ``wavelength`` / ``intensity`` (and common aliases) are preferred.
   Headerless two-column files are also accepted: the leftmost column is
   treated as wavelength and the next as intensity, and a warning is shown.

3. Once the data file is selected a dialogue will open allowing the user to
   specify the intensity units of the spectral data so the GUI can convert
   into the required units for the model. The wavelengths __must__ be in units of nm.

4. Enter a name for the spectrum. It is saved under ``data/user_library/spectra/``
   and appears in the **Spectrum Library**.

The same import dialogue is also available from **Import Spectrum** in the
**Spectrum Library** window.

### Import Stimulus

**Import Stimulus** loads a saved stimulus protocol and registers it in the stimulus
library so it can be selected for a model run.

1. Select **File** then **Import Stimulus**.

2. Choose a ``.json`` stimulus file.

3. If a stimulus with the same name already exists, choose a different name
   to import a copy, or cancel.

4. The stimulus is written to ``data/user_library/stimuli/`` and added to the
   run dropdown, and selected on the **Main Page**.

## Configure

This menu allows the user to modify all of the internal model parameters through
**Model parameters...**. The exact nature of the parameters is detailed in the
Appendix of [`Melanopsin_Model_Tutorial.ipynb`](Melanopsin_Model_Tutorial.ipynb).

Changed model parameters are not saved by default: the user must
click the **Update** button at the bottom right of the pop up **Model parameters** window 
to enable the new parameters on the next model run.

Custom model parameters can be saved by clicking the **Save Parameters** button in the bottom
left hand corner of the **Model parameters** window. Saved custom parameter
sets can also be imported into the GUI by utilizing the **Load Preset** button
that is also in the bottom left hand corner of the **Model parameters** window.

Configuration presets are saved as ``.json`` files and can be found in
``data/user_library/configurations``. That folder includes two shipped
manuscript presets: ``GS.json`` (gain-shifting kinetics, no feedback or
bleaching) and ``No-GS.json`` (the Emanuel and Do 2015 R-M-E tristable model with no
adaptation). Load either with **Load Preset**, then click **Update**. These
are distinct from **Reset to defaults**, which restores the GUI's defaults: the augmented gain shifting (aGS) model with feedback and bleaching on.

To reset the parameters back to their default manuscript values, click the **Reset to defaults** button at the bottom right of the **Model parameters** window.

The model parameters automatically reset to their default values
upon re-launch of the GUI.

## Data

The **Data** menu provides access to saved light spectra and light stimuli. It also allows users to create their own spectra and stimuli.

### Stimulus Builder

The **Stimulus Builder** allows the user to combine multiple built-in and custom spectra
to create a custom stimulus to be run in the model.

The **Stimulus Builder** window is broken up into three sections. The top half of the page is
the **Stimulus Declaration**. This is where stimulus intervals are declared — time intervals in which light has a specific, unchanging spectrum and intensity. Continuous spectral or intensity changes can be approximated by combining consecutive short intervals with small changes. Intervals can also consist of darkness. 

In the bottom left of the **Stimulus Builder** window is the **Light Monitor**.
The **Light Monitor** shows the light intensity (in log10) of the stimulus during each stimulus interval. The currently selected interval is highlighted on the monitor. Clicking an interval on the monitor selects that column in the stimulus table. Left and Right arrow keys move between intervals when you are not editing a text field.

In the bottom right of the **Stimulus Builder** is the **Spectral Viewer**.
The **Spectral Viewer** displays a plot of the light spectrum within a stimulus interval. The user can select which interval is in view of the **spectral viewer** by hovering or clicking an interval
in the **Light Monitor**, selecting an interval in the **Stimulus Declaration**, or using the Left/Right arrow keys.

To create a stimulus interval

1. Click the box in the **Spectrum** row that says **--none--**. This will
   open a dialogue to have the user select a spectrum from the **Spectrum Library**.

2. Click the box in the **Log10 Intensity** row and input a value to represent the total photon intensity in     photons/µm²/s (log₁₀ by default). To express intensity in linear units (instead of log), uncheck the **Log10 intensity input** checkbox at the top-center of the **Stimulus Builder** window.

3. Click the box in the **Duration** row and input a value in seconds. This will be
   the duration that the chosen **Spectrum** will be displayed at the selected **Intensity**.

This is the process for creating a single interval. To create additional intervals
click the **Add...** button in the top left corner to append an additional
blank interval to the end of the stimulus. Right-click an interval number in the
table and choose **Insert before** or **Insert after** to insert a blank interval
at that position. To delete an interval select the
**Delete...** button directly next to **Add...**. This will open a dialogue for
the user to delete a single specific interval.

Once a stimulus has been designed to the user's liking the **Save** button, located
directly next to the **Delete...** button, can be used to save the custom stimulus
and add it to the stimulus library. The ``.json`` of a saved stimulus can be found
in the ``data/user_library/stimuli`` folder in the program files.

Additionally, existing stimuli can be imported into the **Stimulus Builder** for editing. 

- **Load from library...** opens a dialogue allowing the user to select and
   open a stimulus from the **Stimulus Library**. If there are intervals already
   declared in the **Stimulus Declaration** the user will be prompted to select
   whether they want to replace or append the selected stimuli.

- **Load manuscript...** produces a similar dialogue to **Load from library...**
   but only allows the user to select from the manuscript stimuli. This is to
   reduce clutter.

Lastly, the total **Duration** of the stimulus is determined by the sum
of each interval's duration. However, the user has the ability to extend the
**Duration** of the stimulus. If the custom **Duration** is longer than the sum of the intervals then the stimulus is padded with a **dark** interval so that the sum of the intervals matches the custom **Duration**. The **dark**
interval is a built-in spectrum that has zero intensity for all wavelengths.

### Stimulus Library

The stimulus library contains all of the saved stimuli currently in the GUI.
Selecting a stimulus on the left of the window will
allow that stimulus to be viewable on the right of the
window. A search box above the stimulus list filters it by name.
If the stimulus has multiple spectral intervals the user will be
given a digital scrubber to view the respective spectral intervals and
their relative intensities.

To import a new stimulus the user can select the **Import Stimulus**
button in the bottom left hand corner of the window. The same import dialogue 
from **Import Stimulus** in the **File** menu will appear.

To create a new stimulus from scratch, click **Make Custom Stimulus**
in the bottom left hand corner of the **Stimulus Library** window.
This opens the stimulus builder used elsewhere in the GUI.

The stimulus library also allows the user to remove stimuli from the program.
To delete a stimulus

1. Select a stimulus from the stimulus list. Doing so will highlight the stimulus

2. Click **Delete stimulus** in the bottom left hand corner
   of the **Stimulus Library** window.

3. The user will be prompted if they actually want to delete the stimulus from the program.

Deleted stimuli cannot be recovered and will have to be remade/reimported.

The user cannot delete manuscript stimuli from the program. These stimuli
are hard-coded into the program and will always be available upon relaunch
of the GUI.

### Spectrum Library

The spectrum library functions very closely to the stimulus library.
Here the user can select spectra from the spectrum list to view on the right side of the window. A search box above the spectrum list filters it by name.
From this window the user can view, add, delete, and crop spectra.

To add a spectrum the user can

1. Click the **Import Spectrum** button in the bottom left hand corner of the
   **Spectrum Library** window.

2. This will open the same dialogue as described in **Import Spectrum**.

To delete a spectrum the user can

1. Select the spectrum they wish to delete. Doing so will highlight the spectrum.

2. Click the **Delete Spectrum** button in the bottom left hand corner of the
   **Spectrum Library** window.

3. If no saved stimulus protocols reference the spectrum, confirm the deletion
   in the yes/no dialog.

4. If one or more stimulus protocols use the spectrum, a warning dialog lists
   those protocols and asks how to proceed before the spectrum is removed:

   - **Delete these stimulus protocols** — remove every affected protocol from
     the stimulus library, then delete the spectrum.
   - **Replace it in these protocols with:** — pick another spectrum from the
     library; every interval that referenced the doomed spectrum is rewritten
     to the replacement (choosing **Dark** also forces intensity to zero), then
     the original spectrum is deleted.

Deleted spectra cannot be recovered and will have to be remade/reimported.

The user cannot delete manuscript spectra from the program. These spectra
are hard-coded into the program and will always be available upon relaunch
of the GUI.

The user also has the ability to crop spectral data. This is useful when
the extreme ends of measured spectra (very short or very long wavelengths) are noisy.

1. Select the spectrum to crop. Doing so will highlight the spectrum.

2. Once the target spectrum is selected click the **Crop Spectrum** button
   at the bottom of the window next to **Delete Spectrum**. This will
   open up the **Crop Spectrum** window.

3. The user has two methods for cropping the spectrum. First, the user
   can manually select beginning/end wavelengths at the bottom of the window
   by typing in their values.

4. Instead of manually inputting wavelengths to define the interval the user
   can instead click and drag the two vertical bars on the left and right of the
   spectrum. The blue-highlighted spectral data between these two bars will be
   retained while data outside this range will be removed.

5. Once the new spectral interval has been selected the user can save this as
   a new spectrum by clicking the **Save as new spectrum...** button in the bottom
   left hand side of the **Crop Spectrum** window.

6. The user will be prompted to name the new spectrum. Once named, the new spectrum
   will automatically be added to the **Spectrum Library**.

### Spectrum Builder

The **Spectrum Builder** allows the user to create a
**Monochromatic Spectrum** and mix spectra to create a **Mixed Spectrum**.

To make a custom **Monochromatic Spectrum**

1. Select the **Monochromatic Spectrum** tab in the top left hand corner of
   the **Spectrum Builder** window.

2. In the **Target Wavelength** text box input the desired wavelength. The
   **Monochromatic impulse** graph updates automatically as you type.

3. Save the custom spectrum to the **Spectrum Library** by clicking the
   **Save to library...** button in the bottom right hand corner of
   the **Spectrum Builder** window.

The saved **Monochromatic Spectrum** will be saved in the directory
``data/user_library/spectra``.

If the user wants to, instead, create a more complex spectrum they can
select the **Mixed Spectrum** tab directly next to the **Monochromatic Spectrum**
tab at the top of the **Spectrum Builder** window. In this window the user can
add spectra from the **Spectrum Library** to the current custom
spectrum by clicking the **Add row** button at the top of the **Mixed Spectrum** tab.

If the user wishes to delete one of these added rows they can simply
utilize the **Remove last row** button at the top of the **Mixed Spectrum** tab.

After adding a new row (and thus equivalently adding a new light source to the spectrum)
the user can set the **weight** of individual sources using the **weight** text box.

After adding the appropriate weights, the user can choose how each row is scaled
before mixing with the checkboxes at the top of the **Mixed Spectrum** tab:

- **Normalize rows to peak** (selected by default) divides each row by its maximum so
  the row peaks at 1.
- **Normalize rows to integrated photon count** divides each row by its trapezoid
  integral over wavelength so the row integrates to 1.

Only one of these checkboxes can be selected at a time. With neither selected, each
row uses its raw spectral values. Row weights are applied after any normalization.
The **Weighted mixture** graph updates automatically when rows, weights, or
normalization options change.

To save the newly-created custom spectrum, click the **Save to library...** button in the bottom right hand side of the **Mixed Spectrum** window. The **Mixed Spectrum** will be saved in the directory ``data/user_library/spectra``.

## Plot

**Compare saved predictions...** opens the **Prediction Comparison** window, where
several model runs can be overlaid on the same axes. Use **Add files...** to load
previously exported runs, **Load current run** to include the run currently shown on
the **Main Page**, and **Remove selected** to drop a run from the comparison.

**Data Comparator...** overlays recorded experimental data on a model run. Datasets
may be `.csv`, `.tsv`, `.txt`, or Excel (`.xlsx`/`.xls`). If a column named like
`time` is present, its values become the data time base (units such as `ms`, `s`,
or `min` are inferred from the header text and converted to seconds; the trace is
zero-based to match the model). When no time column is found, you are prompted for
a sampling rate (Hz) or sample interval (s). Headered files with multiple data
columns offer a **Data column** picker (default: first non-time column). Legacy
headerless single-column recordings still work and use the sampling-rate prompt
(prefilled at 100 Hz). The model run can be the current run or a previously
exported `.csv` or `.npz`. Recorded data are interpolated onto the model's time
base before plotting. The model series shown defaults to `currentGlobalGain` and
can be changed in the dialogue.

## Settings

The **Settings** menu holds **Auto-save options**.

**Auto-save options** turns auto-save on or off. The user can choose to have the
output model run figure and/or the simulation data saved to a chosen location. The default locations are
`outputs/images` for figures and `outputs/predictions` for simulation data;
**Choose folder...** picks a different location and **Reset to default** restores
the original one. When a checkbox is on, the current save path is displayed.

By default both options are turned off as both have the potential of being very large
in data size depending on the duration of the experiment.

Autosave options are reset to default (off) after relaunch of the GUI.

## Help

The **Help** menu opens project documentation:

- **Open README…** opens the local `README.md` in the default app.
- **Open GitHub…** opens the project repository page in the default browser.




