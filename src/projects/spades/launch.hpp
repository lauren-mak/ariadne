//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/config_struct.hpp"

#include "pipeline/graph_pack.hpp"
#include "stages/construction.hpp"
#include "pipeline/genomic_info_filler.hpp"
#include "gap_closer.hpp"
#include "stages/simplification.hpp"
#include "mismatch_correction.hpp"
#include "pair_info_count.hpp"
#include "second_phase_setup.hpp"
#include "repeat_resolving.hpp"
#include "distance_estimation.hpp"
#include "hybrid_aligning.hpp"
#include "chromosome_removal.hpp"
#include "series_analysis.hpp"
#include "barcode_index_construction.hpp"
#include "pipeline/stage.hpp"
#include "contig_output_stage.hpp"
#include "scaffold_graph_construction_stage.hpp"
#include "scaffolder_analysis_stage.hpp"
#include "barcode_deconvolution_stage.hpp"

// (New) Plan A: gmapper
#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/paths/bidirectional_path_io/bidirectional_path_output.hpp"
#include "io/graph/gfa_reader.hpp"
#include "pipeline/graph_pack.hpp"

namespace spades {

inline bool MetaCompatibleLibraries() {
    const auto& libs = cfg::get().ds.reads;
    if (not (libs[0].is_paired()) or libs[0].is_mate_pair())
        return false;
    if (libs.lib_count() > 2)
        return false;
    if (libs.lib_count() == 2 &&
        libs[1].type() != io::LibraryType::TSLReads &&
        libs[1].type() != io::LibraryType::PacBioReads && libs[1].type() != io::LibraryType::NanoporeReads)
            return false;
    return true;
}

inline bool HybridLibrariesPresent() {
    for (size_t lib_id = 0; lib_id < cfg::get().ds.reads.lib_count(); ++lib_id) 
        if (cfg::get().ds.reads[lib_id].is_hybrid_lib()) 
            return true;
    return false;
}

void assemble_genome() {
    INFO("SPAdes started");
    if (cfg::get().mode == debruijn_graph::config::pipeline_type::meta && !MetaCompatibleLibraries()) {
        ERROR("Sorry, current version of metaSPAdes can work either with single library (paired-end only) "
                      "or in paired-end + (TSLR or PacBio or Nanopore) mode.");
        exit(239);
    }

    StageManager SPAdes({cfg::get().developer_mode,
                         cfg::get().load_from,
                         cfg::get().output_saves});

    size_t read_index_cnt = cfg::get().ds.reads.lib_count();

    debruijn_graph::conj_graph_pack conj_gp(cfg::get().K,
                                            cfg::get().tmp_dir,
                                            read_index_cnt,
                                            cfg::get().ds.reference_genome,
                                            cfg::get().flanking_range,
                                            cfg::get().pos.max_mapping_gap,
                                            cfg::get().pos.max_gap_diff);
    conj_gp.kmer_mapper.Attach();

    if ( cfg::get().run_assembly ) {
        SPAdes.add<debruijn_graph::Construction>()
              .add<debruijn_graph::GenomicInfoFiller>();
        SPAdes.add<debruijn_graph::ContigOutput>();
    } else {
        gfa::GFAReader gfa(cfg::get().assembly_graph);
        gfa.to_graph(conj_gp.g, true);
    }

    SPAdes.add<debruijn_graph::BarcodeDeconvolutionStage>();

    INFO("Starting from stage: " << cfg::get().entry_point);     
    SPAdes.run(conj_gp, cfg::get().entry_point.c_str());

    INFO("SPAdes finished");
}

}
