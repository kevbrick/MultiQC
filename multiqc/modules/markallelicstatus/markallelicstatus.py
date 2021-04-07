#!/usr/bin/env python

""" MultiQC module to parse stats output from markallelicstatus """

import logging
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, linegraph, heatmap
from multiqc.utils import report
from multiqc import config

import os
import re
import numpy as np
from copy import copy
from itertools import combinations_with_replacement, zip_longest
from operator import itemgetter

from .utils import read_allelicstatus_stats, allelictypes_common

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """This MultiQC module parses various
    stats produced by the modified markAllelicStatus.py from HiC Pro."""

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Mark_allelic_status', anchor='markallelicstatus',
        href="https://raw.githubusercontent.com/kevbrick/HiC-Pro/master/scripts/markAllelicStatus.py",
        info="Marks reads mapping to one or other parental genome"
            " Modified from HiCPro to output percentages by type")

        # Find and load any markallelicstatus stats summary files:
        self.markallelicstatus_stats = dict()
        for f in self.find_log_files('markallelicstatus', filehandles=True):
            s_name = f['s_name']
            self.markallelicstatus_stats[s_name] = self.parse_markallelicstatus_stats(f)

        # Filter to strip out ignored sample names
        self.markallelicstatus_stats = self.ignore_samples(self.markallelicstatus_stats)

        if len(self.markallelicstatus_stats) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.markallelicstatus_stats)))

        # # Add to self.js to be included in template
        # self.js = {'assets/js/multiqc_pairtools.js' : os.path.join(os.path.dirname(__file__), 'assets', 'js', 'multiqc_pairtools.js') }
        #
        # # determine max total reads for general stats:
        # self.max_total_reads = 0
        # for s_name in self.markallelicstatus_stats:
        #     self.max_total_reads = \
        #         max(self.max_total_reads, self.markallelicstatus_stats[s_name]['total'])
        # # self.max_total_reads = max(self.pairtools_stats[s_name]['total'] for s_name in self.pairtools_stats)

        #self.markallelicstatus_general_stats()

        # Report sections
        self.add_section (
            name = 'Allelic alignment status',
            anchor = 'pair-types',
            description="Number of read classified according to their parental alignment status",
            helptext = '''For further details check out HiCPro
                        <a href=\"https://github.com/nservant/HiC-Pro\" > HiC Pro</a>
                        documentation.''',
            plot = self.markallelicstatus_chart()
        )

    def parse_markallelicstatus_stats(self, f):
        """ Parse an allelicstatus summary stats file """
        # s_name = f['s_name']
        # f_name = f['fn']
        # log.info("parsing {} {} ...".format(s_name,f_name))
        f_handle = f['f']
        return read_allelicstatus_stats(f_handle)


    def markallelicstatus_chart(self):
        """ Generate the allelic status report """

        report_field = "allelic_status"

        # Construct a data structure for the plot: sample ->
        # and keep track of the common pair types - for nice visuals:
        atypes_dict = dict()
        rare_atypes = set()
        for s_name in self.markallelicstatus_stats:
            atypes_dict[s_name] = dict()
            for atype in self.markallelicstatus_stats[s_name][report_field]:
                atypes_dict[s_name][atype] = self.markallelicstatus_stats[s_name][report_field][atype]
                # # update the collection of rare pair types :
                # if atype not in allelictypes_common:
                #     aare_ptypes.add(atype)

        # organize pair types in a predefined order, common first, rare after:
        atypes_annotated = OrderedDict()
        for atype, color in allelictypes_common.items():
            atypes_annotated[atype] = {'color': color, 'name': atype}

        print (atypes_dict)
        # Config for the plot
        config = {
            'id': 'allelic_status',
            'title': 'mark allelic status: allelic status report',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }

        return bargraph.plot(atypes_dict, atypes_annotated, pconfig=config)
