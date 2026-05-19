# NeumuacFQ — CFTR Precision Medicine Platform

**Version:** 5.1 — updated 19 May 2026  
**Status:** Active development · Target launch June 2026  
**Architecture:** Single-file HTML application · Supabase backend · Deno Edge Functions

---

## What This Is

NeumuacFQ is a clinical-grade CFTR variant intelligence platform built for CF centres. It allows clinicians to look up CFTR variants, understand ETI eligibility, validate classifications, curate missing data, and generate clinical PDF reports — all from a single browser-based interface with no installation required.

The platform is designed around two distinct user modes:

- **Clinical mode** — clean lookup interface for pulmonologists, CF nurses, and geneticists consulting on patient cases
- **Curator mode** — full data management interface for curators actively filling classification gaps and reviewing AI-generated evidence

---

## Technology Stack

| Layer | Technology |
|-------|------------|
| Frontend | Single HTML file — vanilla JS, CSS custom properties, IBM Plex Sans/Mono, Syne |
| Database | Supabase (PostgreSQL) |
| Auth/Write protection | Supabase RPC PIN verification |
| AI evidence scanning | Supabase Edge Functions (Deno) + NCBI eutils |
| PDF generation | Browser `window.print()` — no external library |
| Fonts | Google Fonts (IBM Plex Sans, IBM Plex Mono, Syne) |
| Icons | Font Awesome 6 |

---

## Database Schema

### Core tables

| Table | Purpose |
|-------|---------|
| `variants` | Master variant registry — 2,200+ CFTR variants |
| `evidence_links` | Evidence chain per variant — PubMed, CFTR2, functional studies |
| `validation_history` | Clinical validation records with clinician identity and confidence |
| `variant_notes` | Persistent clinical notes thread per variant |
| `institutions` | CF centres — name, city, country, coordinates |

### Intelligence tables

| Table | Purpose |
|-------|---------|
| `search_logs` | Every search — hit/miss, session, location, institution |
| `search_miss_queue` | Variants searched but not found — human flagging + robotic scanning |
| `search_history_log` | Cross-session search history per institution |
| `helios_pending_review` | Curator review queue — AI proposals + human miss flags |
| `helios_audit_log` | Full audit trail of every Helios decision |
| `helios_bot_runs` | Robotic scan run history |
| `helios_monitor_targets` | Variants under active surveillance |

### Classification tables

| Table | Purpose |
|-------|---------|
| `cftr_class_rules` | Class-level clinical rules (ETI response, modulator eligibility) |
| `contradictions` | Auto-detected data contradictions per variant |
| `cftr_variants_archive` | Classification version history — every change recorded |
| `audit_issues` | Field-level audit flags |
| `research_gaps` | Gap score analysis results |

### Reference tables

| Table | Purpose |
|-------|---------|
| `journal_impacts` | Impact factors for evidence quality scoring |
| `author_publications` | Author h-index cache |
| `evidence_sources` | Trusted source registry |
| `validator_roles` | Clinician role registry |
| `app_config` | Runtime configuration key-value store |

---

## New in Version 5.1 (19 May 2026)

### Forms & data capture
- **Validation modal** — confidence level selector (Confirmed / High / Moderate / Low), institution auto-filled from session, variant state snapshot at time of validation, session ID linkage
- **Evidence modal** — DOI field, publication year, patient population size, session and institution attribution
- **Manual import form** — classification source dropdown (CFTR2, ClinVar, Internal Curation, etc.), attributed to field, clinical phenotype notes
- **Session tracking** — `STATE.sessionViews` accumulates every variant viewed with timestamp, entry method, and app mode

### Export system
- **Export Centre** — full modal replacing the old floating dropdown. Four categories: Clinical Documents, Data & Curation Exports, Research & Governance, Session Report
- **Pre-download identity modal** — name, role, and report purpose captured before every PDF. Persisted to localStorage. Creates audit trail linking PDF to session
- **Elite PDF generator** — IBM Plex Sans/Mono + Syne typography, dark header band, hero section with class block and confidence strip, reporter identity line, DOI and population size in evidence chain, phenotype notes in recommendation card
- **Session Clinical Report** — new export type. Documents every variant viewed, validated, and searched in a session. Suitable for patient encounter notes or audit trails

### Miss system
- **Human miss card** — replaces silent robotic failure. Textarea (500 chars), auto-focused on miss, saves note + institution + session to `search_miss_queue`, creates `HUMAN_MISS_FLAG` item in Helios immediately
- **Robotic lie killed** — fake "searching robotically" notification removed. CORS probe added. Status only set to `in_review` when findings actually confirmed
- **Supabase Edge Function** — server-side PubMed + ClinVar scanning on miss insert. No CORS restriction. Creates `MISS_RESOLVED` Helios item with real paper titles, journals, years, URLs

### UI & modes
- **Clinical / Curator mode toggle** — header segment pill. Persisted to localStorage. Clinical mode hides: Helios button, Analytics button, write-mode lock, Miss/Gaps/Done chips, gap score bar, missing pip indicators, curation sprint hero. Curator mode shows everything
- **Hospital selector redesign** — dropdown with sticky search bar, country grouping, query highlighting, verified badge, usage count, type icon per institution
- **Export Centre** — session context bar shows institution, variant count, views this session, mode, session start time
- **Logo bug fixed** — stray `>` character and `mix-blend-mode:multiply` removed from hospital modal logo

### Detail panel
- **Validation history** — confidence badge (colour-coded), institution with hospital icon, variant state snapshot showing class/ETI/determination at moment of validation
- **Evidence chain** — DOI in blue monospace, population size in green with people icon, both pulled from new DB columns
- **Classification section** — classification source, attributed to, phenotype notes, provenance strip (institution, source, creation date)

### Helios inbox
- **HUMAN_MISS_FLAG** — registered as a type with amber flag icon and "CLINICIAN" badge. Structured reasoning parsed: clinician note shown prominently in italic, institution with hospital icon, session ID in monospace
- **Richer greeting** — decomposes pending count by type: urgent, clinician flags, new evidence
- **Card headers** — relative timestamps, larger type icons with border ring, confidence bar 5px with label
- **Item body** — dark background wash, 16px padding

### Annotation modal
- New editable fields: `classification_source` (dropdown), `attributed_to` (text), `phenotype_notes` (textarea), `alt_names`, `source_summary`
- `cftr_class` now includes Class VI and II/III
- `class_subtype` now includes `presumed` and `atypical`
- `classification_source` and `class_confidence` trigger the reason-required flow

### Curation sprint
- Queue expanded from 5 to 10 variants
- Card redesigned: gap score badge (colour-coded by severity), context strip (searches, validations, confidence, provenance), missing field bullets with hover effects
- Button handlers moved from inline onclick to ID-based wiring (reliable in closures)

### Mobile viewport
- New breakpoints at 768px and 480px
- Hospital dropdown full-width, mode toggle compact, export centre full-screen, compare tray stacks vertically, validation confidence grid 2×2 on small screens, Helios full-screen

### Write-mode PIN
- `id="writeModeConfirmBtn"` added (was missing — caused silent failures)
- Label in separate `<span>` so icon is preserved during loading state
- Slide-in animation on every open (not just first)
- Backdrop blur increased

### Compare tray
- Redesigned header: icon, subtitle, verdict badge with colour dot
- Gradient background
- Allele column labels with dot indicator
- Reads `class_confidence` for display

### Bug fixes
- `Delta Helix CFTR Precision Medicine Platform` replaced with `NeumuacFQ` everywhere
- `deltahelix_hospital` localStorage key preserved (changing it would wipe saved hospitals)
- `cftr_variants_archive` query fixed from `?id=eq.` to `?variant_id=eq.`
- All 404/401/400 Supabase errors diagnosed and resolved via schema migration

---

## Schema Migration — v5.1

Run in Supabase SQL Editor. Safe to run multiple times (all guards use `IF NOT EXISTS`).

See `supabase-edge-function-miss-scanner.md` for the Edge Function deployment and webhook configuration.

Key additions:
- `validation_history` — 8 new columns including `confidence_level`, `institution_name`, `snapshot_cftr_class`, `snapshot_eti`
- `evidence_links` — `doi`, `population_size`, `session_id`, `added_by_institution`
- `variants` — `classification_source`, `attributed_to`, `phenotype_notes`, `added_by_institution`, `import_session_id`
- `search_miss_queue` — `human_flagged`, `flagged_by`, `flagged_at`, `robotic_attempted`, `robotic_blocked`
- `helios_pending_review` — `HUMAN_MISS_FLAG` added to `item_type` CHECK constraint
- `helios_audit_log` — `session_id`, `app_mode`, `institution_name`, `institution_id`
- New tables: `search_history_log`, `variant_notes`
- Fixed: `contradictions` RLS, `cftr_variants_archive` `variant_id` FK

---

## Edge Function — Miss Queue Scanner

**Function name:** `scan-miss-queue`  
**Trigger:** Database webhook on INSERT to `search_miss_queue`  
**Behaviour:** Server-side PubMed + ClinVar scan. Creates `MISS_RESOLVED` Helios item only when findings confirmed. Sets `status: 'in_review'` only then. Records `robotic_attempted` and `robotic_blocked` honestly regardless of outcome.

Full deployment instructions in `supabase-edge-function-miss-scanner.md`.

---

## File Structure

```
index.html                              — entire application (single file)
README.md                               — this file
supabase-edge-function-miss-scanner.md — Edge Function code + deployment guide
```

---

## Key Configuration

Inside `index.html`, find the `CONFIG` object near the top of the JavaScript section:

```javascript
const CONFIG = {
    SUPABASE: {
        URL:      'https://your-project.supabase.co',
        KEY:      'your-anon-key',
        ANON_KEY: 'your-anon-key',
    },
    WRITE_CODE:     'your-pin',        // Fallback PIN if RPC not deployed
    WRITE_DURATION: 30 * 60 * 1000    // Write mode duration (30 min)
};
```

The PIN is verified server-side via `check_write_pin` Supabase RPC. The `WRITE_CODE` is a client-side fallback only — used if the RPC hasn't been deployed.

---

## Known Remaining Work

| Item | Priority | Notes |
|------|----------|-------|
| Full testing session | High | First session with all v5.1 changes live against real DB |
| Write-mode PIN UX | Medium | Gate position and timeout UX — functional but not polished |
| Compare tray | Medium | Functional, redesigned header — content could be richer |
| Supabase Edge Function testing | High | Deployed but not yet tested with a real miss insert |
| `check_write_pin` RPC | Medium | Server-side PIN verification — `migration_002.sql` not yet deployed |
| Light mode | Low | Exists but not fully tested with all new components |

---

## Platform Identity

**Name:** NeumuacFQ  
**Full name:** NeumuacFQ CFTR Precision Medicine Platform  
**Institution:** Neumología · CHUAC  
**localStorage namespace:** `helix_*`, `neumuacfq_*`  
**Session storage key:** `neumuacfq_filter`  
**Legacy storage key:** `deltahelix_hospital` — preserved for backward compatibility, not user-visible
