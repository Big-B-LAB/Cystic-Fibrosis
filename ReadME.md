# NeumuacFQ — CFTR Precision Medicine Platform
## Technical Documentation v4.0

---

## Table of Contents

1. [System Overview](#1-system-overview)
2. [Architecture](#2-architecture)
3. [Data Model](#3-data-model)
4. [Core Algorithms](#4-core-algorithms)
5. [Search Engine](#5-search-engine)
6. [Helios Evidence Engine](#6-helios-evidence-engine)
7. [Robotic Autonomous Search](#7-robotic-autonomous-search)
8. [Compound Heterozygosity Engine](#8-compound-heterozygosity-engine)
9. [ACMG/AMP Classification Framework](#9-acmgamp-classification-framework)
10. [Clinical Reference System](#10-clinical-reference-system)
11. [Clinical Reports & MDT Output](#11-clinical-reports--mdt-output)
12. [Search Intelligence Layer](#12-search-intelligence-layer)
13. [Geographic Intelligence](#13-geographic-intelligence)
14. [Security Model](#14-security-model)
15. [Performance Characteristics](#15-performance-characteristics)
16. [Deployment](#16-deployment)
17. [Roadmap](#17-roadmap)

---

## 1. System Overview

NeumuacFQ is a clinical-grade CFTR variant curation and precision medicine decision support platform for Cystic Fibrosis (CF). It maintains a curated database of CFTR gene variants and provides clinicians with real-time access to classification data, ETI (Elexacaftor/Tezacaftor/Ivacaftor) response predictions, compound genotype analysis, ACMG/AMP pathogenicity scoring, evidence chains, and autonomous AI-driven evidence surveillance.

The platform serves three primary clinical functions:

**Variant lookup** — clinicians search for specific variants by any nomenclature (legacy name, HGVS protein notation, HGVS cDNA notation, or alternative names) and receive immediate classification, clinical reference data, modulator eligibility, and ACMG scoring.

**Compound genotype analysis** — two alleles entered, full clinical prediction output including phenotype severity, sweat chloride range, pancreatic sufficiency probability, per-modulator eligibility, newborn screening detection, and one-click MDT report generation.

**Variant curation** — authorised users curate and enrich variant records. All curation is subject to automated consistency checking by the Helios Evidence Engine. The robotic search engine autonomously queries PubMed, ClinVar, and CFTR2 when a variant is not found, bringing literature directly to the curator's inbox.

### Key Statistics (Production)
- **2,237** CFTR variants under active curation
- **1,339** complete records (class, ETI prediction, HGVS notation)
- **18,700+** lines of code · **359** functions · **1.3 MB** self-contained
- **16** active database tables
- Live clinical use across multiple CF centres in Spain and internationally
- Real-time geographic search intelligence across all sessions

---

## 2. Architecture

### 2.1 Technology Stack

The platform is a single-page application (SPA) delivered as a static HTML file. All persistent state is managed through a Supabase backend.

| Layer | Technology |
|---|---|
| Frontend | Vanilla JavaScript (ES2020+), HTML5, CSS3 |
| Charts | Chart.js (CDN) |
| Geographic map | Leaflet.js + MarkerCluster (CDN) |
| Backend | Supabase (PostgreSQL + PostgREST REST API) |
| Hosting | GitHub Pages or any static file server |
| Offline cache | IndexedDB (browser-native) |
| Fonts | IBM Plex Sans, IBM Plex Mono, Syne (Google Fonts) |
| External APIs | NCBI E-utilities (PubMed + ClinVar), ipapi.co (geolocation) |

### 2.2 Application Architecture

```
CONFIG              — Runtime configuration (Supabase URL, UI constants, write PIN)
STATE               — Global application state (variants, search index, view state)
Database Layer      — loadVariants(), CRUD via Supabase REST, dbWrite() wrapper
Search Engine       — buildSearchIndex(), searchVariants(), normaliseQueryWithSteps()
Nomenclature Layer  — AA3TO1 map, IVS/rsID/chromosomal normalisation, alias table
View Layer          — renderVariantList(), renderVariantDetail()
Clinical Layer      — renderClinicalReferencePanel(), renderACMGPanel()
Compound Het Engine — computeCompoundHetVerdict(), MODULATOR_MAP
Genotype Calculator — openCompoundGenotypeCalculator(), generateMDTReport()
Intelligence Layer  — Search logging, gap scoring, geographic analysis
Helios Engine       — Autonomous consistency checking, evidence surveillance
Robotic Search      — triggerRoboticSearchOnMiss(), multi-source concurrent scan
Offline Engine      — HELIX_CACHE (IndexedDB), failed-save queue
```

### 2.3 State Management

All application state is held in a single `STATE` object that exists only in memory for the duration of a session:

```javascript
STATE = {
    variants: [],               // Full variant array loaded at init
    variantsById: Map,          // O(1) lookup by numeric ID
    variantsByLegacyName: Map,  // O(1) lookup by legacy name string
    searchIndex: {
        lookup: Map,            // normalised term → Set<variantId>
        stats: {}
    },
    masterStats: {},            // Pre-computed aggregate statistics
    view: {
        filteredVariants: [],   // Current filtered/searched result set
        selectedVariant: null,  // Currently displayed variant
        searchQuery: '',        // Live search query
        currentPage: 1,
        sortAsc: true,
        dataViewMode: 'all',
        classFilter: 'all'
    },
    searchIntelligence: {
        sessionId: string,      // sess_ + 8 random alphanumeric chars
        locationCache: {},      // IP geolocation (fetched once per session)
        hospital: string,       // Clinician institution (localStorage)
        recentSearches: [],     // Session hit/miss history (max 20)
        globalHits: Map,        // Query → hit count this session
        globalMisses: Map,      // Query → { count, sessions }
        _lastSearchTier: number // 0=exact 1=prefix 2=substr 3=fuzzy 4=phonetic
    },
    writeMode: {
        active: boolean,
        expiry: timestamp,
        timerInterval: handle
    }
}
```

### 2.4 Initialisation Sequence

1. Render loading screen with animated SVG rings + NeumuacFQ logo
2. `loadClassRules()` + `loadJournalImpacts()` — parallel, non-blocking DB-driven config
3. `loadVariants()` — paginated fetch from Supabase (1,000 records per page with evidence links joined)
4. `buildSearchIndex()` — construct in-memory inverted index
5. `calculateMasterStats()` — compute aggregate counts by class, ETI, validation status
6. `initSearch()` — attach input event listeners with dual debounce (180ms UI, 1200ms log)
7. `loadSearchIntelligence()` — fetch historical search logs and miss queue for dashboard
8. `initHospital()` — retrieve or prompt for institution identity from localStorage
9. `fetchLocationOnce()` — single IP geolocation call for geographic search logging
10. `runSurveillanceAlerts()` — silent background check of monitored variants (3s delay)
11. Loading screen fade-out, application ready

Typical initialisation: **1.2–2.8 seconds** depending on variant count and network latency. On Supabase failure, IndexedDB offline cache is served automatically with a clear banner.

---

## 3. Data Model

### 3.1 Core Tables

**`variants`** — Primary clinical record for each CFTR variant.

| Field | Type | Description |
|---|---|---|
| `legacy_name` | text UNIQUE | Traditional mutation name (e.g. F508del, G551D) |
| `protein_name` | text | HGVS protein notation (e.g. p.Phe508del) |
| `cdna_name` | text | HGVS cDNA notation (e.g. c.1521_1523delCTT) |
| `alt_names` | text | Space/comma/semicolon-separated alternative names |
| `cftr_class` | text | Functional class I–VI, II/III compound |
| `class_subtype` | text | true, presumed, exceptional, atypical, standard |
| `eti_prediction` | text | responsive, non_responsive, unknown |
| `final_determination` | text | CF-causing, VUS, non CF-causing |
| `clinical_alert` | text | **Computed column** — contextual treatment guidance |
| `validation_count` | integer | Number of clinician validations |
| `search_hit_count` | integer | Total successful searches for this variant |
| `changed_prev` | text | Whether classification changed from prior version |

The `clinical_alert` field is a PostgreSQL computed column generating contextual treatment warnings (e.g. *"⚠️ EXCEPTIONAL Class I — Consider ETI despite typical Class I rules"*) based on class/subtype combinations. Surfaced as a colour-coded banner in every variant detail view.

**`evidence_links`** — Evidence sources with classification claims.

Each link carries `cftr_class_claim`, `eti_claim`, `quality_score`, `study_type`, `journal`, `year`, `authors`, and `fingerprint`. The claim fields are the primary input to the Helios contradiction engine.

**`search_logs`** — Every search event (hit or miss) with full geographic context.

Captures `query`, `is_hit`, `matched_field`, `hospital_name`, `country`, `city`, `latitude`, `longitude`, `session_id`. Primary data source for all geographic intelligence features.

**`search_miss_queue`** — Aggregated failed search queries.

Maintains `miss_count`, `unique_sessions`, `top_hospital`, and `notes` per query. Updated via PostgreSQL `ON CONFLICT` upsert. Drives robotic search triggering and miss queue intelligence.

### 3.2 Helios-specific Tables

**`helios_pending_review`** — Primary Helios inbox. All AI findings queue here before any DB modification. Item types: `NEW_EVIDENCE`, `CONTRADICTION`, `DRIFT_ALERT`, `MISS_RESOLVED`, `PATTERN_SIGNAL`, `REGIONAL_SIGNAL`, `INTERNAL_LOGIC`, `INTERNAL_NAMING`, `INTERNAL_ORPHAN`, `INTERNAL_DUPLICATE`, `EXTERNAL_GUIDELINE`.

**`helios_audit_log`** — Append-only record of every clinician decision. Stores `review_snapshot` (full JSONB state at decision time), field-level change diff (`db_change_field`, `db_change_old`, `db_change_new`), clinician name and role, timestamp.

**`helios_bot_runs`** — Operational metadata: scan type, status, variants scanned, findings created, error count, duration.

**`helios_monitor_targets`** — Variants under continuous background literature surveillance.

### 3.3 Research & Governance Tables

| Table | Purpose |
|---|---|
| `contradictions` | Detected evidence contradictions, CRITICAL/HIGH/MEDIUM/LOW severity |
| `cftr_variants_archive` | Version snapshot on every classification-critical field change |
| `cftr_class_rules` | Admin-editable class descriptions — DB-driven, not hardcoded |
| `research_gaps` | Computed gap scores (0–100) per variant |
| `validation_history` | Per-variant clinician validation records with role and timestamp |
| `variant_evidence` | Structured evidence linked to `evidence_sources` reliability registry |
| `journal_impacts` | Impact factor cache for evidence quality scoring |
| `institutions` | Verified institution registry for multi-centre attribution |
| `variant_notes` | Persistent clinical notes per variant, visible to all clinicians |

---

## 4. Core Algorithms

### 4.1 Search Index Construction

`buildSearchIndex()` constructs an inverted index (Map: normalised term → Set of variant IDs) at startup. For each variant, indexed keys include:

- Legacy name (lowercased, HGVS prefix stripped)
- Protein name (lowercased, `p.` stripped)
- cDNA name (lowercased, `c.` stripped; nucleotide change portion also indexed separately)
- All alternative names (split on space/comma/semicolon)
- Positional numbers ≥3 digits (supports "508", "551", "1282" numeric searches)

Result: ~15,000–20,000 index entries for a 2,237-variant database. Search time: **< 5ms** in typical browsers.

### 4.2 Gap Scoring

Each variant receives a gap score (0–100) via `calculateGapScore()`:

| Criterion | Weight |
|---|---|
| No CFTR class assigned | 25 |
| No ETI prediction | 20 |
| No final determination | 15 |
| No evidence links | 15 |
| Fewer than 3 evidence sources | 8 |
| No validation + high search demand (≥5 hits) | 12 |
| High search demand (≥10 hits) | 5 |

Gap levels: CRITICAL (≥75) · SEVERE (≥55) · MODERATE (≥35) · MINOR (≥15) · ADEQUATE (<15). Weights configurable via `cftr_class_rules` table without code changes.

### 4.3 Master Statistics

`calculateMasterStats()` computes aggregate counts used throughout the UI — total by class, ETI distribution, exceptional count, validation coverage, data completeness — recomputed after every write operation.

### 4.4 Evidence Quality Scoring

Papers scored 1–10 by `scoreEvidenceQuality()`:

- Base: 5
- Journal impact factor from `journal_impacts` cache: `min(10, round(2 + IF/3))`
- Recency bonus: +1 if published within last 2 years

---

## 5. Search Engine

### 5.1 Nomenclature Intelligence

`normaliseQueryWithSteps()` transforms clinical typing patterns into searchable keys and records each transformation step for UI display (clinician sees: *"searching as f508del — three-letter AA → single letter"*):

**All 20 amino acids, mid-string:**

```
Phe→F  Gly→G  Arg→R  Trp→W  Asn→N  Asp→D  His→H  Lys→K
Leu→L  Pro→P  Thr→T  Val→V  Ile→I  Ser→S  Ala→A  Cys→C
Gln→Q  Glu→E  Met→M  Tyr→Y
```

**Termination codon synonyms** → `x`: Ter, Stop, Opal, Amber, Ochre

**Special handling:**
- Greek/Unicode: `δ`, `Δ` → `del`
- HGVS prefixes stripped: `p.`, `c.`, `r.`, `m.`, `n.`
- IVS notation: `IVS` → `ivs` (preserved for index lookup)
- Chromosomal positions: `chr7:117,559,590` → stripped entirely
- rsID: passed through unchanged (e.g. `rs75527207`)
- Hyphens and spaces removed within variant names

**Canonical alias table (12+):**

| Input | Normalised to |
|---|---|
| `ΔF508`, `deltaf508`, `delta-f508` | `f508del` |
| `phef508del`, `phe508del`, `df508` | `f508del` |
| `W1282Stop`, `W1282Ter` | `w1282x` |
| `G542Stop`, `G542Ter` | `g542x` |
| `R553Stop`, `R1162Stop` | `r553x`, `r1162x` |

### 5.2 Search Tiers

Results sorted by match quality (0 = best):

| Tier | Type | Example |
|---|---|---|
| 0 | Exact match | `G542X` → G542X |
| 1 | Prefix / positional number | `G551` → G551D; `508` → all 508-position variants |
| 2 | Substring | `del` → all deletion variants |
| 3 | Fuzzy Levenshtein ≤ 2 edits | `G55D` → G551D |
| 4 | Phonetic (consonant skeleton) | near-miss typographic variants |

Fuzzy matching only runs when no direct match exists.

### 5.3 Auto-Select Behaviour

**Single result:** immediately selected, detail panel opens, dropdown closes.

**Tier-0 exact match in multiple results:** immediately selected regardless of total count. List jumps to the correct page (e.g. G542X on page 2 of 31 results), row scrolls into view. No manual pagination required.

```javascript
// Condition: top result is an exact match
tier === 0 && topResult.legacy_name.toLowerCase() === normaliseQuery(q)
```

### 5.4 Tier-Aware Search Feedback

Visual feedback fires **only on tier-0 exact match**. Substring, fuzzy, and phonetic matches are visually silent — they are helpful but not confirmation signals.

| Visual signal | Condition |
|---|---|
| **Green sustained edge** on row | Tier-0, data complete and/or validated |
| **Amber sustained edge** on row | Tier-0, data incomplete or unvalidated |
| Header sweep animation | Tier-0 exact match |
| Detail panel entry glow | Tier-0 exact match |
| Amber-to-green gradient sweep | Exceptional variant |
| Full-screen amber wash | Exceptional variant only |
| No signal | Tier 1–4, or browse selection |

Edges persist until the search is cleared. Data quality assessed live:

```javascript
function variantDataQuality(variant) {
    // confirmed → green edge | incomplete → amber edge
    const complete  = hasClass && hasETI;
    const evidenced = evidence_links.length > 0 || validated;
    if (complete && evidenced) return 'confirmed';
    if (complete)              return 'complete';
    return 'incomplete';
}
```

### 5.5 Live Miss Card

When a variant is not found, the search dropdown transforms into a live miss card:

- Robotic scan status: pulsing indicator while searching, confirmation once queued
- Cross-institutional demand: *"3 other sessions searched this term"*
- Similar variants already in the database with edit distance shown
- External links: CFTR2, ClinVar, PubMed, LOVD
- Clinical note field: annotates `search_miss_queue.notes` for curators before they leave

### 5.6 Cross-Session History

Recent successful searches persisted in `localStorage` (max 50 entries). Shown on empty input focus with relative timestamps (*"3h ago"*, *"yesterday"*). One-click re-run.

---

## 6. Helios Evidence Engine

Helios operates on the governance principle: **the AI proposes, the clinician decides.** The bot reads freely from the database but writes only to `helios_pending_review`. All modifications to `variants` require explicit clinician approval.

### 6.1 Inbox Tabs

| Tab | Contents |
|---|---|
| Inbox | All pending review items, filterable by urgency / item type |
| Scanner | Manual scan controls: internal scan, external PubMed, single variant |
| Validation queue | Variants awaiting clinical validation |
| Calibration | Bot confidence feedback and accuracy tracking over time |
| Audit | Complete action log with field-level diffs |
| Health check | Database integrity diagnostics |

### 6.2 Item Types Generated

| Type | Trigger |
|---|---|
| `NEW_EVIDENCE` | PubMed papers found for unclassified or poorly-evidenced variants |
| `CONTRADICTION` | Evidence `cftr_class_claim` or `eti_claim` contradicts current classification |
| `DRIFT_ALERT` | Classification may be outdated based on new literature |
| `MISS_RESOLVED` | Robotic scan found literature for a searched-but-missing variant |
| `REGIONAL_SIGNAL` | Variant searched from 3+ institutions, unvalidated or unclassified |
| `INTERNAL_LOGIC` | Class I marked ETI responsive without exceptional documentation |
| `INTERNAL_NAMING` | HGVS notation format violations (missing `p.` or `c.` prefix) |
| `INTERNAL_ORPHAN` | No CFTR class, searched 3+ times by clinicians |
| `EXTERNAL_GUIDELINE` | New guideline publication detected |

### 6.3 Evidence Contradiction Detection

`heliosCheckEvidenceContradictions()` reads all `evidence_links` with populated `cftr_class_claim` or `eti_claim`. For each variant:

1. Group claims by variant
2. Compare against current `cftr_class` / `eti_prediction`
3. Weight by source `quality_score`
4. Contradicting weight > confirming weight → priority `URGENT`
5. **Persists to `contradictions` table** (not just in-memory) at CRITICAL/HIGH/MEDIUM severity

### 6.4 Internal Consistency Scan

`heliosCheckVariantConsistency()` evaluates every variant against logical rules — Class I + ETI responsive without exceptional marking, missing HGVS prefixes, unclassified variants with high search demand, duplicate name detection.

### 6.5 Literature Surveillance

`runSurveillanceAlerts()` runs 3 seconds after page load:

1. Fetches `helios_monitor_targets` not checked in 7 days
2. PubMed search per monitored variant (400ms between requests, respecting rate limits)
3. Filters for papers since `last_checked_at`
4. New papers → `NEW_EVIDENCE` Helios item
5. Updates `last_checked_at`

### 6.6 Clinician Decision Flow

```
Helios finding created (URGENT or ROUTINE)
       ↓
Clinician opens inbox → reviews: bot reasoning, confidence %, proposed change vs current
       ↓
       ├── Approve           → applies DB change, writes audit entry
       ├── Approve with note  → same, requires written clinical justification
       ├── Mark Exceptional   → writes class_subtype='exceptional',
       │                        suppresses future flags for this variant
       ├── Reject             → status='rejected', note recorded
       └── Defer              → status='deferred', deferred_until=tomorrow
                                defer_count incremented
                                auto-resurfaced when date passes
```

### 6.7 Calibration System

Reads `helios_audit_log` approved/rejected decisions. Computes per-`item_type` precision rates. Adjusts confidence scoring for future findings. Cached per session.

---

## 7. Robotic Autonomous Search

When a clinician searches for a variant not in the database, NeumuacFQ automatically triggers a multi-source robotic scan — **no curator button required**.

### 7.1 Full Flow

```
Clinician types unknown variant name
       ↓
persistSearchEvent() → miss logged to search_miss_queue (ON CONFLICT upsert)
       ↓ (500ms delay, completely non-blocking)
triggerRoboticSearchOnMiss(normQuery)
       ↓
Three sources scanned concurrently via Promise.allSettled():
   roboticPubMedSearch()   — 8 results via NCBI esearch + esummary
   roboticClinVarSearch()  — 5 clinical records via NCBI clinvar
   roboticCFTR2Search()    — 3 CFTR2-reference publications via targeted PubMed
       ↓
scoreEvidenceQuality()      — weight papers by journal_impacts table
deriveClassificationProposal() — class hint from title keyword analysis
       ↓
INSERT into helios_pending_review — pre-filled dossier with all evidence attached
       ↓
Clinician notification: "NeumuacFQ is searching PubMed · ClinVar · CFTR2 robotically.
                         Check Helios inbox within 24h."
       ↓
Curator one-click Approve → variant created in DB, miss queue resolved,
                             version archived in cftr_variants_archive
```

**Time from miss to live variant: under 24 hours with no manual literature search.**

### 7.2 Rate Limiting

NCBI Entrez: 3 req/s unauthenticated. Concurrent (not sequential) scanning minimises wall time. Each request uses `AbortSignal.timeout(8000)`. Duplicate scan prevention via `ROBOTIC_SCAN_QUEUE` Set.

---

## 8. Compound Heterozygosity Engine

Accessible via the **Genotype** button in the main header toolbar.

### 8.1 Modulator Eligibility Map

```javascript
MODULATOR_MAP = {
    'Trikafta (ETI)':    ≥1 F508del allele OR both alleles ETI-responsive,
    'Kalydeco (IVA)':    ≥1 Class III or IV allele,
    'Symdeko (TEZ/IVA)': ≥1 F508del allele,
    'Orkambi (LUM/IVA)': BOTH alleles F508del (homozygous only)
}
```

23 variants with FDA/EMA ivacaftor monotherapy approval: G551D, G1244E, G1349D, G178R, G551S, S1251N, S1255P, S549N, S549R, R117H + 13 additional gating mutations.

### 8.2 Dominant Allele

Class rank: I=6, II=5, III=4, IV=3, V=2, VI=1. Higher rank = worse function = dominant phenotype driver.

### 8.3 Calculator Output

For two variant names (full database autocomplete):

- Phenotype severity (classic CF · moderate-severe · mild/CFTR-RD)
- Expected sweat chloride range (from dominant allele class)
- Residual CFTR function estimate
- Pancreatic sufficiency probability
- Age of onset
- CFTR-related disorder association
- Newborn screening detection across all three panels
- Per-modulator eligibility for all four approved modulators
- ETI mismatch warning when alleles conflict
- One-click **Generate MDT Report**

If either allele is not in the database, robotic search fires automatically.

### 8.4 Compare Tray

Pin any variant → select a second → compound het tray appears at the bottom of screen with the clinical verdict and per-modulator eligibility panel. Updates live as the second variant changes.

---

## 9. ACMG/AMP Classification Framework

Automatically computed for every variant. Based on Richards et al., *Genetics in Medicine* 17:405–423 (2015).

### 9.1 Criteria

**Pathogenic evidence (auto-assessed where data available):**

| Code | Criterion | Auto-trigger condition |
|---|---|---|
| PVS1 | Null variant, LOF disease mechanism | Class I + CF-causing determination |
| PS3 | Functional study: damaging effect | evidence_link with study_type='functional' |
| PS4 | Prevalence elevated in affected individuals | validation_count ≥ 3 |
| PM2 | Absent/very low frequency in population | population_frequency < 0.001 |
| PM4 | In-frame indel | name contains del/ins/dup, non-Class I |
| PP3 | Computational evidence: damaging | Classes I, II, or III |
| PP5 | Reputable source reports pathogenic | source_summary includes 'CFTR2' |

**Benign evidence:**

| Code | Criterion | Auto-trigger condition |
|---|---|---|
| BA1 | Allele frequency > 5% | population_frequency > 0.05 |
| BS1 | Allele frequency > 1% | population_frequency > 0.01 |
| BP4 | Computational: no impact | ETI non-responsive + no class |

### 9.2 Classification Thresholds

| Score | Classification |
|---|---|
| ≥ 10 | Pathogenic |
| 6–9 | Likely pathogenic |
| 0–5 | Uncertain significance |
| −6 to −1 | Likely benign |
| ≤ −7 | Benign |

Displayed in every variant detail with met criteria highlighted. Score and classification also appear in the Clinical PDF and MDT Report. Clearly labelled as automated — expert review required before clinical reporting.

---

## 10. Clinical Reference System

Every variant detail panel shows class-specific reference data from hardcoded tables (no additional DB queries):

### 10.1 Class Reference Data

| Class | Sweat Cl⁻ | Residual CFTR | Pancreatic sufficiency | Onset | Primary modulator |
|---|---|---|---|---|---|
| I | ≥80 mmol/L | <1% | <5% | Neonatal | ETI if ≥1 responsive allele |
| II | ≥80 mmol/L | 0–3% | <10% | Neonatal/infancy | ETI strongly indicated |
| III | 60–100 mmol/L | 1–5% | 10–25% | Infancy–childhood | Ivacaftor / ETI |
| IV | 30–60 mmol/L | 5–30% | 30–70% | Childhood–adulthood | Ivacaftor may benefit |
| V | 30–70 mmol/L | 3–10% | 40–80% | Childhood–adulthood | ETI uncertain |
| VI | 40–80 mmol/L | 5–25% | 25–60% | Variable | Under investigation |

All labelled as population-level estimates, not patient-specific.

### 10.2 Newborn Screening Panels

| Panel | Variants |
|---|---|
| UK CLAPA | 27 variants: F508del, G542X, G551D, R117H, N1303K, W1282X, R553X + 20 others |
| ACMG 23+4 | 23 ACMG core variants |
| EU Consensus | 17 pan-European variants |

Panel membership shown per variant: *"This variant is included in X national NBS panel(s) — would be detected at birth."* Non-panel variants explicitly noted.

### 10.3 Ivacaftor vs ETI Distinction

The panel explicitly distinguishes ivacaftor monotherapy eligibility from Trikafta eligibility — these are different clinical indications frequently conflated in practice. *"Ivacaftor monotherapy approved"* shown for the 23 FDA/EMA-labelled gating mutations, separate from the Class II ETI indication.

---

## 11. Clinical Reports & MDT Output

### 11.1 Clinical Variant PDF

Generated by `exportClinicalReport()` — single variant or batch (all classified). Opens in a print-optimised popup window (A4, correct margins).

**Contents:** NeumuacFQ logo · variant name block with badges · classification details + ACMG score · clinical reference grid (sweat chloride, residual function, pancreatic sufficiency, onset, CFTR-RD, modulator) · NBS panel inclusion · clinical recommendation · validation status · PubMed references · legal footer.

### 11.2 MDT Report

Generated by `generateMDTReport()` from the Compound Genotype Calculator. Structured as a clinical letter for MDT correspondence.

**Contents:** NeumuacFQ logo · institution + date · genotype section (both alleles with ACMG classification) · clinical prediction grid (6 fields) · therapeutic eligibility table (per modulator + Trikafta eligibility + dominant allele) · ACMG/AMP table for both alleles · NBS detection status · legal disclaimer (*"For MDT discussion only. Not a substitute for clinical genetics assessment or specialist genetics review."*).

### 11.3 Research Exports

**Filtered cohort CSV** — exports current filter state, not the full database. Includes gap scores and clinical alerts. For grant applications and literature reviews.

**IRB data package (JSON)** — full provenance per variant: nomenclature, classification, therapeutic, validation, evidence, gap analysis, demand signal. For IRB submissions.

**Search Intelligence PDF** — auto-generated research report with hit rate analysis, miss queue, top variants, geographic distribution, institutional context.

---

## 12. Search Intelligence Layer

### 12.1 Session Architecture

Each session receives a random `session_id` (`sess_` + 8 alphanumeric characters) — provides session-level aggregation without tracking individuals or requiring authentication.

### 12.2 Dual Debounce

```
Keystroke → 180ms  → UI filter updates (always responsive)
Keystroke → 1200ms → search_logs INSERT + search_miss_queue UPSERT
```

Prevents every partial keystroke from generating a log entry while keeping the UI instantaneous.

### 12.3 Geographic Enrichment

On first search, `fetchLocationOnce()` calls `ipapi.co/json/` once — city-level geolocation, cached for the full session. All subsequent log entries in the session carry the same geographic context at zero additional API cost.

### 12.4 Dashboard Tabs

**Overview** — Hit rate, unique queries, field breakdown (which nomenclature type is searched most), session sparkline, coverage bar, miss queue chip.

**Miss queue** — All unresolved misses sorted by cross-institutional demand. Misses from 3+ institutions auto-escalate to Helios URGENT.

**Geography** — Leaflet map with marker clusters, choropleth overlay (hit rate by country), period filters (7d / 30d / all time), per-location hit rate bars.

**Velocity heatmap** — Hour-of-day × day-of-week search volume from `search_logs`. Green = high hit rate. Blue = high miss rate. Shows when clinical demand peaks — useful for staffing and surveillance scheduling.

**Co-search clusters** — Variant pairs searched together in the same session. Identifies potential compound het pairs from real clinical behaviour. Cross-institutional pairs highlighted. Click any pair to load both into the compare tray.

**Gap queue** — Variants ranked by research gap score with configurable weights.

---

## 13. Geographic Intelligence

### 13.1 Data Structure

Every `search_logs` entry: `country`, `region`, `city`, `latitude`, `longitude`, `hospital_name`. Supports analysis from individual hospital level to country level.

### 13.2 Cross-Institutional Demand Signal

`heliosRunGeoSignal()` identifies variants searched from multiple institutions in recent logs. A variant searched from 5 different countries warrants priority validation and evidence enrichment even if it's already classified — this is clinically active data. For HIGH signals (3+ hospitals or 3+ countries) on unvalidated/unclassified variants, a `REGIONAL_SIGNAL URGENT` Helios item is created automatically.

### 13.3 Co-Search Intelligence

`renderCoSearchClusters()` reads session co-occurrences — variants searched together in the same session. Groups by session, builds pair counts, ranks by frequency, identifies compound het candidates from real usage patterns. Directly actionable: click any pair to load both into the compound het compare tray.

---

## 14. Security Model

### 14.1 Supabase Row Level Security

RLS enabled on all tables. The `sb_publishable_*` key in the HTML grants read access to variant data and write access to search logs and Helios queues. It does not grant direct write access to `variants`.

All variant modifications go through Helios's review approval flow or through `saveAnnotation()` which requires write-mode PIN activation. This is an application-layer constraint — for production with a broader user base, a backend proxy pattern should replace the client-side key.

### 14.2 Write Mode

PIN-gated, 30-minute sessions. All writes go through `dbWrite()`:
- 3-retry with exponential backoff (800ms × attempt)
- Error classification: network / validation / permission / duplicate / server
- Failed-save queue in `localStorage` — queued writes replay automatically on reconnect
- Connection indicator in header (green / amber-saving / red-error dot)

### 14.3 Audit Trail

**`helios_audit_log`** — append-only. Every decision stores: `review_snapshot` (full JSONB state), `db_change_field` / `db_change_old` / `db_change_new`, `clinician_name`, `clinician_role`, `created_at`.

**`cftr_variants_archive`** — version snapshot triggered on every change to `cftr_class`, `eti_prediction`, `final_determination`, `class_subtype`, or `eti_evidence_level`. Full classification history is recoverable at any point in time.

### 14.4 Formal Dispute Workflow

The dispute button (⚖) on every variant allows a second clinician to formally challenge a classification. Requires: name/role, proposed correct value, evidence/reasoning. Writes to `contradictions` at CRITICAL severity and escalates to Helios URGENT immediately. Audit logged.

---

## 15. Performance Characteristics

### 15.1 Variant Loading

Paginated batches of 1,000 with evidence links joined inline (`select=*,evidence_links(*)`). A 2,237-variant database: ~3 requests at 300–600ms each on Supabase free tier.

### 15.2 Search Performance

O(1) Map lookups for the inverted index. ~15,000–20,000 index entries for 2,237 variants. Complete search (all tiers, including phonetic): **< 5ms** in typical browsers. Levenshtein fuzzy matching only runs when no direct match is found.

### 15.3 Helios Scan Performance

Internal consistency scan yields every 20 variants (`await setTimeout(50ms)`) to avoid blocking the event loop. Full 2,237-variant scan: **5–8 seconds** at this throttle.

External PubMed scans: 350ms between requests (NCBI Entrez limit: 3 req/s unauthenticated; register an API key for 10 req/s).

### 15.4 Calibration Load

Two parallel requests on Helios open: `helios_audit_log` (limit 500) + `helios_pending_review` (limit 500). Both resolve in < 400ms typically. Cached per session.

### 15.5 Offline Mode

Full dataset cached to IndexedDB on every successful load. On startup failure: serves from cache with offline banner. Write operations queue locally in `localStorage` key `helix_failed_saves_v1` and replay on reconnect.

---

## 16. Deployment

### 16.1 Current Setup

- Single `neumuacFQ.html` (~18,700 lines, 1.3 MB including embedded logo as base64 JPEG)
- Hosting: GitHub Pages or any static HTTPS file server
- HTTPS required (ipapi.co geolocation requires secure context)

### 16.2 Configuration

```javascript
// neumuacFQ.html — CONFIG block
const CONFIG = {
    SUPABASE: {
        URL:  'https://[project-ref].supabase.co',
        KEY:  'sb_publishable_[key]'    // publishable key with RLS enabled
    },
    UI: {
        ITEMS_PER_PAGE:  25,
        SEARCH_DELAY:    180,    // ms — UI filter debounce
        LOG_DELAY:       1200,   // ms — search event log debounce
        MIN_SEARCH_CHARS: 2
    },
    WRITE_CODE:     '[PIN]',             // change before deployment
    WRITE_DURATION: 30 * 60 * 1000      // 30 minutes
};
```

### 16.3 Supabase Setup

Required tables (see Section 3): `variants`, `evidence_links`, `search_logs`, `search_miss_queue`, `validation_history`, `helios_pending_review`, `helios_audit_log`, `helios_bot_runs`, `helios_monitor_targets`, `research_gaps`, `contradictions`, `cftr_variants_archive`, `cftr_class_rules`, `institutions`, `journal_impacts`, `variant_notes`

CORS origins — Supabase Dashboard → API Settings:
```
https://[your-domain]
```

### 16.4 Known Limitations

**Single file architecture** — 18,700+ lines. Works well at current scale. Will benefit from ES module decomposition (`Vite` or native `import`) as the team grows beyond one developer.

**Client-side key** — publishable key visible in HTML source. Acceptable with RLS for a trusted clinical team; replace with a backend proxy for any public-facing deployment.

**Session identity** — institution name in `localStorage` is self-reported, not authenticated. For regulated clinical use requiring individual audit trails, integration with an identity provider (NHS login, hospital SSO, SAML 2.0 / OIDC) is required.

**PubMed rate limits** — 3 unauthenticated req/s. Register an NCBI API key for 10 req/s if higher-volume robotic scanning is needed.

---

## 17. Roadmap

Three features that would make NeumuacFQ citable in a peer-reviewed methods paper:

**1. Functional study data layer**
Structured fields for Ussing chamber measurements, organoid swelling assays, Fischer rat thyroid cell assays. Currently the ACMG PS3 criterion (functional study) is inferred from `evidence_link.study_type` — it should be first-class structured data with quantitative values and assay metadata.

**2. Population stratification**
R117H with 5T/7T/9T poly-T tract and TG repeat context has radically different clinical penetrance — currently modelled as a single variant. Haplotype modifier fields and ancestry-stratified allele frequencies from gnomAD are required to model this correctly.

**3. Classification history UI**
`cftr_variants_archive` writes a snapshot on every classification change. A timeline view in the variant detail panel showing *when* a variant moved from VUS to CF-causing, *what evidence* triggered it, and *which Helios approval* authorised the change would provide the audit trail needed for clinical governance and peer review.

---

## Acknowledgements

Built on CFTR2 variant database reference data (Johns Hopkins University / University of North Carolina). ACMG/AMP classification criteria per Richards et al., *Genetics in Medicine* 17:405–423 (2015). Modulator eligibility per current EMA/FDA prescribing information. NBS panel composition per CLAPA UK, ACMG, and European Cystic Fibrosis Society consensus. Geographic data via ipapi.co. Map tiles via CartoDB.

---

*NeumuacFQ v4.0 · Servicio de Neumología · CHUAC — Complexo Hospitalario Universitario de A Coruña*
*For clinical and research use · Not for patient distribution without clinician sign-off*
*Document version: 2.0 · Platform version: 4.0.0-PRODUCTION · Last updated: April 2026*
