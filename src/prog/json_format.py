# == Native Modules
import textwrap
# == Installed Modules

# == Project Modules


def json_format():
	print('''
=== CUSTOM BATCH INSTRUCTIONS ===
The current message contains instructions to guide the
correct formatting of custom batch files to allow simultaneous
ingestion of several editors/base-editors by "medit guide_prediction --custom-batch".

== Custom Batch Template [Base-Editors] ==
*= Explicit category names =* 
The custom base-editors batch content is based on other parameters described in "medit guide_prediction --help", such as ‘--pam’, ‘--pamisfirst’,’--guidelen’,’--edit_win’, ’--target_base’, and ’--result_base’

{
	"<Base Editor Name>": [
		["--pam", --pamisfirst, --guidelen, [--edit_win], "<keep empty>"],
		["--target_base--result_base", "<Base Editor Name>"]
	]
}

*= Example Json file content: CBE, and ABE =*

{
	"CBE": [
		["NGG", false, 20, [4, 8], ""],
		["CT", "CBE"]
	],
	"ABE": [
		["NGG", false, 20, [4, 8], ""],
		["AG", "ABE"]
	]
}

== Custom Batch Template [Endonucleases] ==
*= Explicit category names =*
The custom endonucleases batch content is based on other parameters described in "medit guide_prediction --help", such as ‘--pam’, ‘--pamisfirst’,’--guidelen’, and '--dsb_pos'

{
	"<Endonucleases Name>": [
		'--pam', --pamisfirst, --guidelen, --dsb_pos, "<keep empty>"	
	]
} 

*= Example Json file content: spCas9, and saCas9 =*

{
	"spCas9": [
		'NGG', False, 20, -3, "<keep empty>"	
	],
	"spCas9": [
		'NNGRRT', False, 21, -3, '<keep empty>'
	]
}	
	''')
