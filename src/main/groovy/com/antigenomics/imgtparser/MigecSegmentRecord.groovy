/*
 * Copyright 2013-2015 Mikhail Shugay (mikhail.shugay@gmail.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.antigenomics.imgtparser

class MigecSegmentRecord {
    final String species, gene, segment, id, sequence
    final int referencePoint, cdr1Start, cdr1End, cdr2Start, cdr2End

    MigecSegmentRecord(String species, String gene, String segment, String id, String sequence) {
        this(species, gene, segment, id, sequence, -1)
    }

    MigecSegmentRecord(String species, String gene, String segment, String id, String sequence,
                       int referencePoint) {
        this(species, gene, segment, id, sequence, referencePoint, -1, -1, -1, -1)
    }

    MigecSegmentRecord(String species, String gene, String segment, String id, String sequence,
                       int referencePoint, int cdr1Start, int cdr1End, int cdr2Start, int cdr2End) {
        this.species = species
        this.gene = gene
        this.segment = segment
        this.id = id
        this.sequence = sequence
        this.referencePoint = referencePoint
        this.cdr1Start = cdr1Start
        this.cdr1End = cdr1End
        this.cdr2Start = cdr2Start
        this.cdr2End = cdr2End
    }

    final static String HEADER = "#species\tgene\tsegment\tid\t" +
            "reference_point\tsequence\t" +
            "cdr1.start\tcdr1.end\tcdr2.start\tcdr2.end"

    @Override
    String toString() {
        [species, gene, segment, id, referencePoint, sequence,
         cdr1Start, cdr1End, cdr2Start, cdr2End].join("\t")
    }
}