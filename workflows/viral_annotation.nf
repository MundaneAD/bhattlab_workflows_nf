/* VIRAL ANNOTATION
 * Runs VIBRANT on contigs, runs CheckV on the viruses,
 * and then dereplicates the viral set across samples.
*/

include { input_assembly } from './modules/input/input_assembly'

include { vibrant } from '../modules/viral_annotation/vibrant'
include { genomad } from '../modules/viral_annotation/genomad'


workflow {
  ch_contigs = input_assembly()
  ch_contigs.view()
  ch_viruses = vibrant(ch_contigs)
}

workflow.onComplete {
	log.info ( workflow.success ? "\nViral annotation is done! Yay!\n" : "Oops .. something went wrong" )
}
