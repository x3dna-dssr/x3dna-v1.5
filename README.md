# 3DNA (v1.5) — Historical Reference Release

[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](LICENSE)
[![Status](https://img.shields.io/badge/Status-Public_Archive-inactive.svg)]()
[![Manual](https://img.shields.io/badge/User_Manual-PDF-red.svg)](x3dna-v1.5.pdf)

**3DNA (v1.5)** (released November 2002) is the foundational version of the **3DNA** suite of programs for the analysis, rebuilding, and visualization of 3D nucleic acid structures, associated with the original 2003 *Nucleic Acids Research* publication.

> **Historical Note & Modern Successors:**  
> This version of 3DNA is outdated and no longer actively supported, though it continues to function as originally designed. It is preserved on GitHub as an open-source archive for historical reference and for third-party software pipelines that specifically depend on the v1.5 layout (such as **3D-DART**). Its base-pair identification algorithm (`find_pair`) also served as the underlying engine for **RNAView**.  
>  
> For the final version of the classic 3DNA suite, see **3DNA v2.4**. For modern nucleic acid structural analysis, modeling, and visualization, **3DNA** has been superseded by [**DSSR**](http://x3dna.org/) and the [**wDSSR** web application](https://web.x3dna-dssr.org).

---

## Key Features & Mathematical Foundation

Written in strict ANSI C with helper Perl scripts, **3DNA v1.5** introduced several foundational capabilities to structural bioinformatics:

* **Standard Base Reference Frame:** Implements the standard base-centered reference frame recommended by the 1999 Tsukuba Workshop and approved by the IUBMB/IUPAC Joint Commission on Biochemical Nomenclature (Olson et al., 2001), using Horn's unit-quaternion least-squares fitting procedure.
* **Rigorous & Reversible Analysis/Rebuilding:** Follows the matrix-based **CEHS** mathematical formalism (El Hassan & Calladine, 1995; Lu et al., 1997) to compute six local base-pair parameters (*Shear, Stretch, Stagger, Buckle, Propeller, Opening*), six dimeric step parameters (*Shift, Slide, Rise, Tilt, Roll, Twist*), and six local helical parameters (*x-displacement, y-displacement, Helical rise, Inclination, Tip, Helical twist*), enabling exact 3D atomic and block reconstruction.
* **Automated Base-Pair & Step Classification:** Automatically identifies canonical and non-canonical base pairs (`find_pair`) and classifies dinucleotide steps as A-, B-, or TA-like based on phosphorus ($Z_p$) positioning.
* **55 Fiber Diffraction Models:** Generates 55 regular DNA and RNA helical structures (`fiber`) for arbitrary sequences following NDB atom naming and ordering conventions.
* **Publication-Quality Schematics:** Produces Calladine–Drew style rectangular base/base-pair blocks (`rebuild`, `alc2img`, `pdb2img`, `blocview`) and standardized base-stacking diagrams (`stack2img`) in PostScript, XFig, Alchemy, and Raster3D formats.

For complete documentation—including a step-by-step mathematical walkthrough of the quaternion least-squares fitting, reference frame construction, and local step/helical rebuilding equations—please refer to the **[3DNA v1.5 User Manual (`x3dna-v1.5.pdf`)](x3dna-v1.5.pdf)**.

---

## Directory Structure & Setup

### Package Organization

* **`BASEPARS/`** — Standard residue geometries (`Atomic_?.pdb` and `ATOMIC/`, `RNA_BASES`), Calladine–Drew block definitions (`Block_*.alc` and `BLOCK/`), phosphate coordinates (`Pxyz.dat`, `PxyzH.dat`), residue mapping (`baselist.dat`), online help (`help3dna.dat`), and plotting parameter files (`misc_3dna.par`, `ps_image.par`, `fig_image.par`, `raster3d.par`, `ndb_raster3d.par`).
* **`FIBER/`** — Repeating unit coordinates (`Str01`–`Str55` and `Data/`) for 55 regular DNA and RNA fiber diffraction models.
* **`Examples/`** — Worked examples covering `Analyze_Rebuild/`, `Calladine_Drew/`, `NMR/`, `Stacking/`, and `Triplex/`.
* **`src/`** — Complete ANSI C source code and `Makefile` for building the 19 compiled executables.
* **`bin/`** — Perl utility scripts (`block_atom`, `blocview`, `cp_std`, `dcmnfile`, `del_ms`, `manalyze`, `mstack2img`, `nmr_ensemble`, `nmr_strs`) and target directory for compiled C executables (`alc2img`, `analyze`, `anyhelix`, `cehs`, `comb_str`, `ex_str`, `fiber`, `find_pair`, `frame_mol`, `get_part`, `o1p_o2p`, `pdb2img`, `r3d_atom`, `rebuild`, `regular_dna`, `rotate_mol`, `stack2img`, `std_base`, `step_hel`).

### Compilation & Environment Configuration

1. **Compile the C executables:**
   ```bash
   cd src
   make
   cd ..
   ```

2. **Configure your shell environment:**  
   Set the environment variable `X3DNA` to point to the installation directory and add `$X3DNA/bin` to your command-line search `PATH`.

   **For `bash` / `zsh` (add to `~/.bashrc` or `~/.zshrc`):**
   ```bash
   export X3DNA=/path/to/x3dna-v1.5
   export PATH=$X3DNA/bin:$PATH
   ```

   **For `csh` / `tcsh` (add to `~/.cshrc`):**
   ```csh
   setenv X3DNA /path/to/x3dna-v1.5
   set path = ($path$X3DNA/bin)
   ```

---

## Citation

If you use or reference **3DNA** in your work, please cite at least one of the primary 3DNA publications:

1. **Primary 3DNA Paper (v1.5):**  
   Lu, X.-J., & Olson, W. K. (2003). 3DNA: a software package for the analysis, rebuilding and visualization of three-dimensional nucleic acid structures. *Nucleic Acids Research*, 31(17), 5108–5121. [https://doi.org/10.1093/nar/gkg680](https://doi.org/10.1093/nar/gkg680)

2. **3DNA Protocols Paper:**  
   Lu, X.-J., & Olson, W. K. (2008). 3DNA: a versatile, integrated software system for the analysis, rebuilding and visualization of three-dimensional nucleic-acid structures. *Nature Protocols*, 3(7), 1213–1227. [https://doi.org/10.1038/nprot.2008.104](https://doi.org/10.1038/nprot.2008.104)

3. **Standard Base Reference Frame (Tsukuba Workshop):**  
   Olson, W. K., Bansal, M., Burley, S. K., Dickerson, R. E., Gerstein, M., Harvey, S. C., Heinemann, U., Lu, X.-J., Neidle, S., Shakked, Z., Sklenar, H., Suzuki, M., Tung, C.-S., Westhof, E., Wolberger, C., & Berman, H. M. (2001). A standard reference frame for the description of nucleic acid base-pair geometry. *Journal of Molecular Biology*, 313(1), 229–237. [https://doi.org/10.1006/jmbi.2001.4987](https://doi.org/10.1006/jmbi.2001.4987)

### BibTeX

```bibtex
@article{Lu2003_3DNA,
  author  = {Lu, Xiang-Jun and Olson, Wilma K.},
  title   = {{3DNA}: a software package for the analysis, rebuilding and visualization of three-dimensional nucleic acid structures},
  journal = {Nucleic Acids Research},
  volume  = {31},
  number  = {17},
  pages   = {5108--5121},
  year    = {2003},
  doi     = {10.1093/nar/gkg680}
}

@article{Lu2008_3DNA_NatProtoc,
  author  = {Lu, Xiang-Jun and Olson, Wilma K.},
  title   = {{3DNA}: a versatile, integrated software system for the analysis, rebuilding and visualization of three-dimensional nucleic-acid structures},
  journal = {Nature Protocols},
  volume  = {3},
  number  = {7},
  pages   = {1213--1227},
  year    = {2008},
  doi     = {10.1038/nprot.2008.104}
}

@article{Olson2001_BaseFrame,
  author  = {Olson, Wilma K. and Bansal, Manju and Burley, Stephen K. and Dickerson, Richard E. and Gerstein, Mark and Harvey, Stephen C. and Heinemann, Udo and Lu, Xiang-Jun and Neidle, Stephen and Shakked, Zippora and Sklenar, Heinz and Suzuki, Masashi and Tung, Chang-Shung and Westhof, Eric and Wolberger, Cynthia and Berman, Helen M.},
  title   = {A standard reference frame for the description of nucleic acid base-pair geometry},
  journal = {Journal of Molecular Biology},
  volume  = {313},
  number  = {1},
  pages   = {229--237},
  year    = {2001},
  doi     = {10.1006/jmbi.2001.4987}
}
```

---

## Community & Support

For questions, comments, or discussions regarding the 3DNA/DSSR ecosystem, please visit [GitHub Discussions](https://github.com/orgs/x3dna-dssr/discussions).

## License

Distributed under the [BSD 3-Clause License](LICENSE).