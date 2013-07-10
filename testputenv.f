	program testputenv
	integer putenv,system

	print *,putenv('AAP=aap')
	print *,putenv('NOOT=noot    ')
	print *,system('printenv')

	end
