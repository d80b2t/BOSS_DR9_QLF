<?php include '../../header.php'; echo head('The Sloan Digital Sky Survey Quasar Catalog: ninth data release'); ?>

<h2 id="intro">The SDSS-DR9 Quasar Catalog (DR9Q)</h2>

<!-- INTRODUCTION -->

<p>
Following the tradition established by SDSS-I/II, the SDSS-III BOSS collaboration produces a quasar catalog.
In addition to the usual information provided by the previous versions of the SDSS Quasar Catalogues
(<em>e.g.</em> <a href="http://adsabs.harvard.edu/abs/2010AJ....139.2360S">Schneider et al. 2010</a>),
this new version contains refined redshifts, Broad Absorption Line (BAL) and also emission line properties.
</p>

<p>
The SDSS-DR9 Quasar catalog is fully described in
<a href="http://arxiv.org/abs/1210.5166">Pâris et al. 2012</a>.
We summarize below the main characteristics of this catalog.
</p>

<!-- CATALOG PRODUCTION -->
<h2 id="production">Catalog production</h2>

<h3>Introduction</h3>

<p>
In order to provide the first measurement of the BAO scale in the Lyman-α forest at z ∼ 2.5, BOSS aims at successfully
observing over 150,000 quasars in the redshift range 2.2 ≤ z ≤ 4, where at least part of the Lyman-α forest lies in the
SDSS spectral range. The target selection has been designed in order to provide at least 15 quasars with z ≥ 2.2 per square degree.
The catalog is therefore not uniform by construction but a uniform sample is also identified.
Since BAO measurements require a catalog of maximal purity, all quasar candidates have been visually inspected.
</p>

<h3>Inspected quasar candidates</h3>

<p>
The objects that were inspected have been selected in three different ways:
</p>
<dl>
<dt>Objects targeted as quasars by the main BOSS survey</dt>
<dd>The quasar target selection is described on the
<a href="dr9/algorithms/boss_quasar_ts.php"> quasar target selection page</a>
and in detail in <a href="http://adsabs.harvard.edu/abs/2012ApJS..199....3R">Ross et al. (2012)</a>.</dd>

<dt>Ancillary programs</dt>
<dd>Some BOSS ancillary programs target quasars through peculiar criteria
(special selection, search for variability, peculiar AGN population etc...).
Those objects have been included in the catalog, and extra flags
indicate which ancillary program they are part of.
The whole description of these flags can be found as part of the
<a href="dr9/algorithms/ancillary/"> DR9 documentation</a> and in
the <a href="http://adsabs.harvard.edu/abs/2012arXiv1208.0022D">BOSS overview paper</a>.</dd>

<dt>Serendipitous objects</dt>
<dd>Quasars might be targeted by chance by other programs, such as the luminous galaxy survey.
In order to be as complete as possible, we also tried to identify serendipitous quasars.
In particular, objects that the pipeline identified as QSO with z&gt;2, or GALAXY/BROADLINE objects have
been inspected.</dd>
</dl>

<h3>Visual inspection</h3>

<p>
The visual inspection is performed to
(i) secure the identification of the objects,
(ii) reliably estimate the emission redshift of quasars, and
(iii) identify peculiar features sucha as Damped Lyman-α systems (DLA) and Broad Absorption Lines (BAL).
We therefore manually confirmed or modified the identification of the object and,
when needed, corrected the redshift provided by the pipeline, <em>i.e.</em> when it was wrong
(when <em>e.g.</em> an emission line is misidentified or a bad feature is
considered an emission line) or inaccurate (when emission lines are correctly
identified but not properly centered).
</p>
<p>
After visual inspection, each quasar candidate is classified among one of the following categories:
QSO, QSO_BAL, QSO_Z? when we know this is a quasar but its redshift is not certain, QSO_? when the object
is possibly a quasar, Star, Star_? when the object is possibly a star, Galaxy, ? when the identification is uncertain or Bad
when the data are not good enough to ascertain identification.
The SDSS-DR9 Quasar Catalog only contains secure identifications (<em>i.e.</em> QSO and QSO_BAL).
A supplemental list of quasars also contains all the objects classified as QSO_Z?.
We also provide in a separate list all the objects classified as QSO_?.
</p>


<!-- CATALOG CONTENT DESCRIPTION -->

<h2>Catalog content</h2>

<p>
The full content of the SDSS-DR9 Quasar Catalog content can be found in Table 4 of <a href="http://arxiv.org/abs/1210.5166">Pâris et al. (2012)</a>
and the detailed description of the data model is available
<a href="http://data.sdss3.org/datamodel/files/BOSS_QSO/DR9Q/DR9Q.html">here</a>.
</p>

<!-- DESCRIPTION OF THE FILES -->

<h2 id="files">Description of the files</h2>

<h3>The SDSS-DR9 Quasar Catalog (DR9Q)</h3>
<p>
This file contains all the quasars of the SDSS-DR9 Quasar Catalog.
87,822 quasars have been identified.
This sample is the one used by the SDSS-BOSS collaboration.
</p>
<p>
<a href="http://data.sdss3.org/sas/dr9/boss/qso/DR9Q/DR9Q.fits">DR9Q.fits</a> (82 MB): SDSS-DR9Q (Main catalog, fits format)
</p>
<p>
<a href="http://data.sdss3.org/sas/dr9/boss/qso/DR9Q/DR9Q.dat">DR9Q.dat</a> (108 MB): SDSS-DR9Q (Main catalog, ascii format)
</p>

<h3>Supplemental list</h3>
<p>
We also provide a supplemental list of an additional 949  quasars that have
been identified, respectively, among the galaxy targets of the BOSS and among
missed QSO targets and identified after freezing the SDSS-DR9 Quasar Catalog.</p>

<p>This file also contains all the quasars classified as QSO_Z? for which the identification is secure but the redshift not accurate.</p>

<p>Though the supplemental list has strictly the same format as
the main catalog, note that pieces of information can be missing for
some objects. In that case, the corresponding value in the column
(<em>e.g.</em> "-1" or "-9999" or "0" etc...) has no meaning (see <a href="http://arxiv.org/abs/1210.5166">Pâris et al. 2012</a>).
Only the identification of the object should be considered as secure.
</p>

<p>
<a href="http://data.sdss3.org/sas/dr9/boss/qso/DR9Q/DR9Q_sup.fits">DR9Q_sup.fits</a> (983 kB): Supplemental list (fits format)
</p>
<p>
<a href="http://data.sdss3.org/sas/dr9/boss/qso/DR9Q/DR9Q_sup.dat">DR9Q_sup.dat</a> (1.3 MB): Supplemental list (ascii format)
</p>


<h3>List of potential quasars (QSO_?)</h3>

<p>
We provide the list of objects identified as QSO_? during the visual inspection.
These objects are likely quasars but their identification can not be certain.
This list contains the first seven columns described in the main catalog only.
</p>

<p>
<a href="http://data.sdss3.org/sas/dr9/boss/qso/DR9Q/DR9Q_unconfirmed.fits">DR9Q_unconfirmed.fits</a> (76 kB): List of QSO_? (fits format)
</p>
<p>
<a href="http://data.sdss3.org/sas/dr9/boss/qso/DR9Q/DR9Q_unconfirmed.dat">DR9Q_unconfirmed.dat</a> (100 kB): List of QSO_? (ascii format)
</p>

<?php echo foot(); ?>

