JULIA?=/Applications/Julia-1.7.app/Contents/Resources/julia/bin/julia

.PHONY: docs

coverage: 
	$(JULIA) -e 'using Pkg; Pkg.activate("."); Pkg.test(;test_args=[""],coverage=true)'
	$(JULIA) -e 'using Coverage; coverage=process_folder(); LCOV.writefile("lcov.info", coverage); clean_folder(".")'
	
clean_mem:
	find . -name "*.[0-9][0-9][0-9][0-9][0-9].mem" -exec rm  {} \;	

docs:
	$(JULIA) --color=yes docs/make.jl	

	
