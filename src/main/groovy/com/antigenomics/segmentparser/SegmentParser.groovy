/*
 * Copyright 2013-2017 Mikhail Shugay (mikhail.shugay@gmail.com)
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


package com.antigenomics.segmentparser

def cli = new CliBuilder(usage: 'SegmentParser [options] input_imgt_genedb output_file_prefix')
cli.b(longOpt: 'report-bad', 'reports \"bad\" IMGT records, i.e. those that are not V/D/J segment or ' +
        'V and J segments that do not have a reference point (conserved Cys or Phy/Trp)')
cli.n(longOpt: 'non-func', 'include non-functional alleles (ORF and Pseudogene) in output')
cli.m(longOpt: 'minor-allele', 'output minor alleles (*02, *03 and so on)')
cli.s(longOpt: 'subspecies', 'more species detalisation (e.g. Balb/c and C57Bl6 for MusMusculus)')
cli.h('display help message')

def opt = cli.parse(args)

if (opt == null || opt.arguments().size() < 2) {
    println "[ERROR] Too few arguments provided"
    cli.usage()
    System.exit(-1)
}

if (opt.h) {
    cli.usage()
    System.exit(0)
}


def imgtInputFileName = opt.arguments()[0], outputFilePrefix = opt.arguments()[1]

if (!new File(imgtInputFileName).exists()) {
    println "[ERROR] Input file not exists"
    System.exit(-1)
}

def reader = new FastaReader(imgtInputFileName)
boolean subspecies = opt.s, nonFunctional = opt.n, minorAlleles = opt.m
def imgtParser = new ImgtToMigecParser(nonFunctional, minorAlleles)

def outputFile = new File(outputFilePrefix + ".txt"),
    outputFileCdr12 = new File(outputFilePrefix + ".cdr12.txt")
if (outputFile.absoluteFile.parentFile)
    outputFile.absoluteFile.parentFile.mkdirs()

outputFile.withPrintWriter { pw ->
    outputFileCdr12.withPrintWriter { pw12 ->
        pw.println(MigecSegmentRecord.HEADER)
        pw12.println("species\tgene\tseqnt\tseqaa\t" +
                "cdr1nt\tcdr2nt\tcdr2.5nt\t" +
                "cdr1aa\tcdr2aa\tcdr2.5aa")
        reader.each { fastaRecord ->
            def imgtRecord = new ImgtRecord(fastaRecord, subspecies)
            def migecRecord = imgtParser.parseRecord(imgtRecord)
            if (migecRecord) {
                pw.println(migecRecord)

                if (migecRecord.segment == "Variable" &&
                        [migecRecord.cdr1Start, migecRecord.cdr2Start,
                         migecRecord.cdr1End, migecRecord.cdr2End].every { it >= 0 }) {
                    def cdr1nt = migecRecord.sequence.substring(migecRecord.cdr1Start, migecRecord.cdr1End),
                        cdr2nt = migecRecord.sequence.substring(migecRecord.cdr2Start, migecRecord.cdr2End),
                        cdr25nt = migecRecord.cdr25Start >= 0 ?
                                migecRecord.sequence.substring(migecRecord.cdr25Start, migecRecord.cdr25End) : "",
                        cdr1aa = Util.translateLinear(cdr1nt),
                        cdr2aa = Util.translateLinear(cdr2nt),
                        cdr25aa = Util.translateLinear(cdr25nt)

                    pw12.println([migecRecord.species, migecRecord.id,
                                  migecRecord.sequence, Util.translateLinear(migecRecord.sequence),
                                  cdr1nt, cdr2nt, cdr25nt,
                                  cdr1aa, cdr2aa, cdr25aa].join("\t"))
                }
            }
        }
    }
}

if (opt.b) {
    new File(outputFilePrefix + ".novrefpoint").withPrintWriter { pw ->
        imgtParser.failedVReferencePoint.each {
            pw.println(it)
        }
    }
    new File(outputFilePrefix + ".nojrefpoint").withPrintWriter { pw ->
        imgtParser.failedJReferencePoint.each {
            pw.println(it)
        }
    }
    new File(outputFilePrefix + ".othersegm").withPrintWriter { pw ->
        imgtParser.otherSegment.each {
            pw.println(it)
        }
    }
}
new File(outputFilePrefix + ".metadata").withPrintWriter { pw ->
    pw.println(ImgtToMigecParser.HEADER)
    pw.println(imgtParser.toString())
}