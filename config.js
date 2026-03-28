// ============================================================
// Δ HELIX v3 — CONFIGURATION
// ============================================================

export const CONFIG = {
    SUPABASE: {
        URL: 'https://oxlkisjwfarseyqgeyuw.supabase.co',
        KEY: 'sb_publishable_vKCe5bc1Qr20_LI_awKIYw_CqSFrblh'
    },
    UI: {
        ITEMS_PER_PAGE: 25,
        SEARCH_DELAY: 180,
        NOTIFICATION_DURATION: 4000,
        MIN_SEARCH_CHARS: 2
    },
    CLASSES: ['I', 'II', 'III', 'IV', 'V', 'VI'],
    CLASS_LABELS: {
        I: 'Nonsense / Frameshift — No functional protein',
        II: 'Misfolded protein — Trafficking defect',
        III: 'Gating defect — Reduced channel opening',
        IV: 'Conductance — Reduced channel conductance',
        V: 'Splicing / Reduced synthesis',
        VI: 'Reduced stability / Accelerated turnover'
    },
    CLASS_COLORS: {
        I: '#c0392b',
        II: '#b7610a',
        III: '#0b6e6d',
        IV: '#1a56a5',
        V: '#6b3fa0',
        VI: '#9b59b6'
    },
    CLASS_BG: {
        I: '#fdf2f2',
        II: '#fef8f0',
        III: '#f2fafa',
        IV: '#f0f5fd',
        V: '#f6f2fd',
        VI: '#f9f0ff'
    }
};

export const TRANSLATIONS = {
    en: {
        // Header
        appTitle: 'CFTR Variant Explorer',
        variants: 'Variants',
        etiResponsive: 'ETI Responsive',
        missingData: 'Missing Data',
        analytics: 'Analytics',
        export: 'Export',
        addVariant: 'Add Variant',
        
        // Sidebar
        searchPlaceholder: 'Search variants, cDNA, protein…',
        all: 'All',
        exceptional: 'Exceptional',
        missing: 'Missing',
        complete: 'Complete',
        variantsCount: 'variants',
        
        // Detail
        class: 'CLASS',
        clinicalClassification: 'Clinical Classification',
        therapeuticPrediction: 'Therapeutic Prediction',
        evidenceChain: 'Evidence Chain',
        clinicalValidation: 'Clinical Validation',
        nomenclature: 'Nomenclature',
        legacyName: 'Legacy Name',
        proteinNotation: 'Protein Notation',
        cdnaNotation: 'cDNA Notation',
        
        // Buttons
        validate: 'Validate',
        edit: 'Edit',
        delete: 'Delete',
        save: 'Save',
        cancel: 'Cancel',
        add: 'Add',
        viewAll: 'View All',
        
        // Empty state
        selectVariant: 'Select a variant',
        emptySub: 'Choose a CFTR variant from the explorer to view comprehensive clinical data, therapeutic predictions, and evidence chain.',
        
        // Form labels
        determination: 'Determination',
        subtype: 'Subtype',
        confidence: 'Confidence',
        evidenceLevel: 'Evidence Level',
        recommendation: 'Recommendation',
        notSpecified: 'Not specified',
        
        // ETI prediction
        etiResponsiveLabel: 'ETI Responsive',
        etiResponsiveDesc: 'Predicted to respond to elexacaftor/tezacaftor/ivacaftor',
        etiNonResponsiveLabel: 'Non-responsive',
        etiNonResponsiveDesc: 'Not predicted to respond to current modulators',
        etiUnknownLabel: 'Unknown',
        etiUnknownDesc: 'Insufficient data to predict ETI response',
        
        // Validation
        yourName: 'Your Name',
        yourRole: 'Your Role',
        notes: 'Notes (optional)',
        validationHistory: 'Validation History',
        totalValidations: 'Total',
        uniqueValidators: 'Unique',
        distinctRoles: 'Roles',
        noValidations: 'Not yet validated by clinical staff',
        firstToValidate: 'Be the first to validate',
        latest: 'LATEST',
        
        // Evidence
        evidenceSources: 'Evidence Sources',
        noEvidence: 'No evidence links yet',
        addEvidence: 'Add Evidence',
        sourceType: 'Source Type',
        title: 'Title',
        url: 'URL',
        sourceName: 'Source Name',
        
        // Modals
        addImportVariant: 'Add / Import Variant',
        jsonImport: 'JSON',
        manualForm: 'Manual',
        bulkCsv: 'Bulk CSV',
        jsonPlaceholder: 'Paste JSON here…',
        validateJson: 'Validate',
        import: 'Import',
        clear: 'Clear',
        saveVariant: 'Save Variant',
        
        // Delete
        deleteVariant: 'Delete Variant',
        deleteWarning: 'This cannot be undone',
        permanentlyDelete: 'Permanently delete',
        andAssociatedData: 'and all associated data?',
        
        // Notifications
        loadingVariants: 'Loading variants…',
        connectingDatabase: 'Connecting to database',
        connectionFailed: 'Connection failed',
        retry: 'Retry',
        saved: 'Saved',
        deleted: 'Deleted',
        exportComplete: 'Export complete',
        
        // Classes descriptions
        classIDesc: 'Nonsense / Frameshift — No functional protein produced',
        classIIDesc: 'Misfolded protein — Trafficking defect',
        classIIIDesc: 'Gating defect — Protein reaches cell surface but fails to open',
        classIVDesc: 'Conductance — Reduced channel conductance',
        classVDesc: 'Splicing / Reduced synthesis',
        classVIDesc: 'Reduced stability / Accelerated turnover'
    },
    
    es: {
        // Header
        appTitle: 'Explorador de Variantes de CFTR',
        variants: 'Variantes',
        etiResponsive: 'Respuesta a ETI',
        missingData: 'Datos pendientes',
        analytics: 'Estadísticas',
        export: 'Exportar',
        addVariant: 'Añadir variante',
        
        // Sidebar
        searchPlaceholder: 'Buscar variantes, cDNA, proteína…',
        all: 'Todas',
        exceptional: 'Excepcional',
        missing: 'Pendientes',
        complete: 'Completas',
        variantsCount: 'variantes',
        
        // Detail
        class: 'CLASE',
        clinicalClassification: 'Clasificación clínica',
        therapeuticPrediction: 'Predicción terapéutica',
        evidenceChain: 'Cadena de evidencia',
        clinicalValidation: 'Validación clínica',
        nomenclature: 'Nomenclatura',
        legacyName: 'Nombre histórico',
        proteinNotation: 'Notación proteica',
        cdnaNotation: 'Notación cDNA',
        
        // Buttons
        validate: 'Validar',
        edit: 'Editar',
        delete: 'Eliminar',
        save: 'Guardar',
        cancel: 'Cancelar',
        add: 'Añadir',
        viewAll: 'Ver todas',
        
        // Empty state
        selectVariant: 'Seleccione una variante',
        emptySub: 'Elija una variante de CFTR del explorador para consultar los datos clínicos, predicciones terapéuticas y cadena de evidencia.',
        
        // Form labels
        determination: 'Determinación',
        subtype: 'Subtipo',
        confidence: 'Confianza',
        evidenceLevel: 'Nivel de evidencia',
        recommendation: 'Recomendación',
        notSpecified: 'No especificado',
        
        // ETI prediction
        etiResponsiveLabel: 'Responde a ETI',
        etiResponsiveDesc: 'Predicción de respuesta a elexacaftor/tezacaftor/ivacaftor',
        etiNonResponsiveLabel: 'No responde',
        etiNonResponsiveDesc: 'No se predice respuesta a moduladores actuales',
        etiUnknownLabel: 'Desconocido',
        etiUnknownDesc: 'Datos insuficientes para predecir respuesta a ETI',
        
        // Validation
        yourName: 'Su nombre',
        yourRole: 'Su cargo',
        notes: 'Notas (opcional)',
        validationHistory: 'Historial de validaciones',
        totalValidations: 'Total',
        uniqueValidators: 'Únicos',
        distinctRoles: 'Cargos',
        noValidations: 'Aún no validada por personal clínico',
        firstToValidate: 'Sea el primero en validar',
        latest: 'MÁS RECIENTE',
        
        // Evidence
        evidenceSources: 'Fuentes de evidencia',
        noEvidence: 'No hay enlaces de evidencia aún',
        addEvidence: 'Añadir evidencia',
        sourceType: 'Tipo de fuente',
        title: 'Título',
        url: 'URL',
        sourceName: 'Nombre de la fuente',
        
        // Modals
        addImportVariant: 'Añadir / Importar variante',
        jsonImport: 'JSON',
        manualForm: 'Formulario',
        bulkCsv: 'CSV múltiple',
        jsonPlaceholder: 'Pegue JSON aquí…',
        validateJson: 'Validar',
        import: 'Importar',
        clear: 'Limpiar',
        saveVariant: 'Guardar variante',
        
        // Delete
        deleteVariant: 'Eliminar variante',
        deleteWarning: 'Esta acción no se puede deshacer',
        permanentlyDelete: 'Eliminar permanentemente',
        andAssociatedData: 'y todos los datos asociados?',
        
        // Notifications
        loadingVariants: 'Cargando variantes…',
        connectingDatabase: 'Conectando con la base de datos',
        connectionFailed: 'Error de conexión',
        retry: 'Reintentar',
        saved: 'Guardado',
        deleted: 'Eliminado',
        exportComplete: 'Exportación completada',
        
        // Classes descriptions
        classIDesc: 'Sin sentido / Cambio de pauta de lectura — Sin proteína funcional',
        classIIDesc: 'Proteína mal plegada — Defecto de tráfico',
        classIIIDesc: 'Defecto de apertura — El canal llega a la superficie pero no se abre',
        classIVDesc: 'Defecto de conductancia — Flujo reducido',
        classVDesc: 'Empalme / Síntesis reducida',
        classVIDesc: 'Estabilidad reducida / Recambio acelerado'
    }
};
