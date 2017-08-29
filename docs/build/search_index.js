var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Welcome!",
    "title": "Welcome!",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Welcome!-1",
    "page": "Welcome!",
    "title": "Welcome!",
    "category": "section",
    "text": "Functionality for working with exact single-arm two-stage designs for clinical trials with binary endpoint.Special focus lies on the computation of optimal two-stage designs extending the ideas of Simon [1] using the Julia programming language. BinaryTwoStageDesigns relies heavily on the functionality of the JuMP package and is tested against the commercial mixed integer programming solver Gurobi. While JuMP supports other solvers the use of Gurobi is highly recommended. Academic licenses for Gurobi can be obtained free of charge upon request.[1] Simon, R. Optimal two-stage designs for phase II clinical trials. Controlled Clinical Trials 1989; 10, 1–10"
},

{
    "location": "index.html#Quickstart-1",
    "page": "Welcome!",
    "title": "Quickstart",
    "category": "section",
    "text": "A quickstart guide is available as jupyter notebook."
},

{
    "location": "index.html#Contents-1",
    "page": "Welcome!",
    "title": "Contents",
    "category": "section",
    "text": "Depth = 3\nPages = [\n  \"sample_space.md\",\n  \"parameters.md\",\n  \"designs.md\",\n  \"optimal_designs.md\"\n]"
},

{
    "location": "index.html#Index-1",
    "page": "Welcome!",
    "title": "Index",
    "category": "section",
    "text": ""
},

{
    "location": "sample_space.html#",
    "page": "Sample Spaces",
    "title": "Sample Spaces",
    "category": "page",
    "text": ""
},

{
    "location": "sample_space.html#BinaryTwoStageDesigns.SampleSpace",
    "page": "Sample Spaces",
    "title": "BinaryTwoStageDesigns.SampleSpace",
    "category": "Type",
    "text": "SampleSpace\n\nDefines restrictions on the sample space of a two-stage design to be used during optimization.\n\n\n\n"
},

{
    "location": "sample_space.html#BinaryTwoStageDesigns.interimsamplesizerange-Tuple{BinaryTwoStageDesigns.SampleSpace}",
    "page": "Sample Spaces",
    "title": "BinaryTwoStageDesigns.interimsamplesizerange",
    "category": "Method",
    "text": "interimsamplesizerange(ss::SampleSpace)\n\nReturn integer vector of valid n_1 values.\n\n\n\n"
},

{
    "location": "sample_space.html#BinaryTwoStageDesigns.maxsamplesize-Tuple{BinaryTwoStageDesigns.SampleSpace}",
    "page": "Sample Spaces",
    "title": "BinaryTwoStageDesigns.maxsamplesize",
    "category": "Method",
    "text": "maxsamplesize(ss::SampleSpace)\n\nReturn maximal overal sample size as integer.\n\n\n\n"
},

{
    "location": "sample_space.html#BinaryTwoStageDesigns.maxsamplesize-Union{Tuple{BinaryTwoStageDesigns.SampleSpace{TI,TR},TI2}, Tuple{TI2}, Tuple{TI}, Tuple{TR}} where TI2<:Integer where TR<:Real where TI<:Integer",
    "page": "Sample Spaces",
    "title": "BinaryTwoStageDesigns.maxsamplesize",
    "category": "Method",
    "text": "maxsamplesize(\n  ss::SampleSpace{TI,TR}, n1::TI2\n) where {TI<:Integer,TR<:Real,TI2<:Integer}\n\nReturn maximal sample size given n1 as integer.\n\n\n\n"
},

{
    "location": "sample_space.html#BinaryTwoStageDesigns.isgroupsequential-Tuple{BinaryTwoStageDesigns.SampleSpace}",
    "page": "Sample Spaces",
    "title": "BinaryTwoStageDesigns.isgroupsequential",
    "category": "Method",
    "text": "isgroupsequential(ss::SampleSpace)\n\nIndicate whether the sample space is restricted to group-sequential designs.\n\n\n\n"
},

{
    "location": "sample_space.html#The-Sample-Space-1",
    "page": "Sample Spaces",
    "title": "The Sample Space",
    "category": "section",
    "text": "Sample space objects can be used to encode the feasible space of binary  two-stage designs for optimization.SampleSpace\n\ninterimsamplesizerange(ss::SampleSpace)\n\nmaxsamplesize(ss::SampleSpace)\n\nmaxsamplesize(ss::SampleSpace{TI,TR}, n1::TI2) where {TI<:Integer,TR<:Real,TI2<:Integer}\n\nisgroupsequential(ss::SampleSpace)"
},

{
    "location": "parameters.html#",
    "page": "Design Parameters",
    "title": "Design Parameters",
    "category": "page",
    "text": ""
},

{
    "location": "parameters.html#Design-Parameters-1",
    "page": "Design Parameters",
    "title": "Design Parameters",
    "category": "section",
    "text": ""
},

{
    "location": "parameters.html#BinaryTwoStageDesigns.Parameters",
    "page": "Design Parameters",
    "title": "BinaryTwoStageDesigns.Parameters",
    "category": "Type",
    "text": "Parameters\n\nAbstract type representing a generic set of paramters for finding optimal two-stage designs. \n\n\n\n"
},

{
    "location": "parameters.html#BinaryTwoStageDesigns.label-Tuple{BinaryTwoStageDesigns.Parameters}",
    "page": "Design Parameters",
    "title": "BinaryTwoStageDesigns.label",
    "category": "Method",
    "text": "label(par::Parameters)\n\nReturn the label of a Parameters-object.\n\n\n\n"
},

{
    "location": "parameters.html#BinaryTwoStageDesigns.null-Tuple{BinaryTwoStageDesigns.Parameters}",
    "page": "Design Parameters",
    "title": "BinaryTwoStageDesigns.null",
    "category": "Method",
    "text": "null(par::Parameters)\n\nReturn the response rate under the null hypothesis of a Parameters-object.\n\n\n\n"
},

{
    "location": "parameters.html#BinaryTwoStageDesigns.mtoer-Tuple{BinaryTwoStageDesigns.Parameters}",
    "page": "Design Parameters",
    "title": "BinaryTwoStageDesigns.mtoer",
    "category": "Method",
    "text": "mtoer(par::Parameters)\n\nReturn the maximal type one error rate of a Parameters-object.\n\n\n\n"
},

{
    "location": "parameters.html#BinaryTwoStageDesigns.mcrv-Tuple{BinaryTwoStageDesigns.Parameters}",
    "page": "Design Parameters",
    "title": "BinaryTwoStageDesigns.mcrv",
    "category": "Method",
    "text": "mcrv(par::Parameters)\n\nReturn the minimal clinically relevant response rate of a Parameters-object.\n\n\n\n"
},

{
    "location": "parameters.html#BinaryTwoStageDesigns.samplespace-Tuple{BinaryTwoStageDesigns.Parameters}",
    "page": "Design Parameters",
    "title": "BinaryTwoStageDesigns.samplespace",
    "category": "Method",
    "text": "samplespace(par::Parameters)\n\nReturn the SampleSpace-object of a Parameters-object.\n\n\n\n"
},

{
    "location": "parameters.html#Overview-1",
    "page": "Design Parameters",
    "title": "Overview",
    "category": "section",
    "text": "Objects of type Parameters are used to define objective functions and their  respective parameters for finding optimal two-stage designs. Every Parameters object also contains a SampleSpace which restricts  the feasible region of the sample space and implements a score method  which is to be optimized.Parameters\n\nlabel(par::Parameters)\n\nnull(par::Parameters)\n\nmtoer(par::Parameters)\n\nmcrv(par::Parameters)\n\nsamplespace(par::Parameters)"
},

{
    "location": "parameters.html#BinaryTwoStageDesigns.MESS",
    "page": "Design Parameters",
    "title": "BinaryTwoStageDesigns.MESS",
    "category": "Type",
    "text": "MESS{T_samplespace<:SampleSpace,TR<:Real} <: PointAlternative\n\nMESS{T_samplespace<:SampleSpace}(\n  samplespace::T_samplespace,\n  p0, p1,\n  alpha, beta,\n  pess;\n  minstoppingforfutility::Real   = 0.0,\n  minconditionalpower::Real      = 0.0,\n  MONOTONECONDITIONALPOWER::Bool = true\n)\n\nThis type represents a set of parameters for finding optimal two-stage designs minimizing the expected sample size on a point in the parameter space subject to type one and two error rate constraints.\n\n> [1] Simon R. Optimal two-stage designs for phase II clinical trials. `Controlled Clinical Trials` 1989; 10, 1-10.\n\n> [2] Kunzmann K and Kieser M. Optimal adaptive two-stage designs for single-arm trials with binary endpoint. `arxive.org` 2016; arXiv:1605.00249.\n\nParameters\n\nParameter Description\nsamplespace a sample space object\np0 upper boundary of the null hypothesis\np1 point alternative to power on\nalpha maximal tolerable type one error rate\nbeta maximal tolerable type two error rate on p1\npess response rate under which to minimize expected sample size\nminstoppingforfutility minimal probability for stopping for futility under p0\nminconditionalpower minimal conditional power upon continuation to stage two\nMONOTONECONDITIONALPOWER if true, the conditional power must be monotonously increasing, this constraint is only relevant if nmax is set very restrictively\n\n\n\n"
},

{
    "location": "parameters.html#MESS-1",
    "page": "Design Parameters",
    "title": "MESS",
    "category": "section",
    "text": "MESS"
},

{
    "location": "parameters.html#BinaryTwoStageDesigns.LiuScore",
    "page": "Design Parameters",
    "title": "BinaryTwoStageDesigns.LiuScore",
    "category": "Type",
    "text": "LiuScore{TI<:Integer,TR<:Real}\n\nLiuScore(\n  samplespace::SampleSpace,\n  p0::Real, pmcrv::Real, prior, \n  alpha::Real, beta::Real,\n  fs::Real, fp::Real;\n  npriorpivots::Integer = 50,\n  npivots::Integer = 15,\n  minconditionalpower::Real = 0.0,\n  MONOTONECONDITIONALPOWER::Bool = false,\n  label::String = \"\"\n)\n\nThis type represents a set of parameters for finding optimal two-stage designs minimizing a variant of the criterion proposed by Liu et al. in [1]. \n\n[1] Liu G, Frank G, Zhu C, and Lu C. Evaluating the adaptive performance of flexible sample size designs with treatment difference in an interval. Statistics in Medicine 2008; 27(4): 584-596.\n\nParameters\n\nParameter Description\nsamplespace a sample space object\np0 upper boundary of the null hypothesis\npmcrv point alternative to power on\nprior prior distribution on true response rate, expected score under the prior conditional on a relevant effect is minimized\nalpha maximal tolerable type one error rate\nbeta maximal tolerable type two error rate on p1\nfs scaling parameter as described in [1]\nfp scaling parameter as described in [1]\nnpriorpivots number of pivot points for integration of prior\nnpivots number of pivot points piecewise linear approximation of local score\nminconditionalpower minimal conditional power upon continuation to stage two\nMONOTONECONDITIONALPOWER if true, the conditional power must be monotonously increasing, this constraint is only relevant if nmax is set very restrictively\nlabel descriptive string for parameter set\n\n\n\n"
},

{
    "location": "parameters.html#LiuScore-1",
    "page": "Design Parameters",
    "title": "LiuScore",
    "category": "section",
    "text": "LiuScore"
},

{
    "location": "parameters.html#BinaryTwoStageDesigns.EB",
    "page": "Design Parameters",
    "title": "BinaryTwoStageDesigns.EB",
    "category": "Type",
    "text": "EB{TI<:Integer,TR<:Real} <: VagueAlternative\n\nEB( # default values\n  samplespace::SampleSpace,\n  p0::Real, pmcrv::Real, prior,\n  gamma::Real, lambda::Real;\n  a::Real = 1, b::Real = 1, targetpower::Real = .8, k::Real = 1,\n  alpha::Real = 0.05,\n  minconditionalpower::Real = 0.0, MONOTONECONDITIONALPOWER::Bool = true,\n  npriorpivots::Integer = 50, ngpivots::Integer = 15,\n  label::String = \"\"\n)\n\nThis type represents a set of parameters for finding optimal two-stage designs maximizing the expected benefit.\n\nParameters\n\nParameter Description\nsamplespace a sample space object\np0 upper boundary of the null hypothesis\npmcrv point alternative to power on\nprior prior distribution on true response rate, expected score under the prior conditional on a relevant effect is minimized\ngamma per patient costs\nlambda benefit of correctly rejecting null if relevant effect exists\na scaling parameter of nonlinear power transfrom (set to infinity for thresholding)\nb scaling parameter of nonlinear power transfrom\ntargetpower target power for thresholding\nk factor for non-responder costs\nalpha maximal type one error rate\nminconditionalpower minimal conditional power upon continuation to stage two\nMONOTONECONDITIONALPOWER if true, the conditional power must be monotonously increasing, this constraint is only relevant if nmax is set very restrictively\nnpriorpivots number of pivot points for integration of prior\nnpivots number of pivot points piecewise linear approximation of local score\nlabel descriptive string for parameter set\n\n\n\n"
},

{
    "location": "parameters.html#EB-1",
    "page": "Design Parameters",
    "title": "EB",
    "category": "section",
    "text": "EB"
},

{
    "location": "parameters.html#BinaryTwoStageDesigns.MBESS",
    "page": "Design Parameters",
    "title": "BinaryTwoStageDesigns.MBESS",
    "category": "Type",
    "text": "MBESS{TI<:Integer,TR<:Real} <: VagueAlternative\n\nMBESS(\n  samplespace::SampleSpace,\n  p0::Real, pmcrv::Real, prior;\n  a::Real = 1, b::Real = 1, targetpower::Real = .8, k::Real = 1,\n  alpha::Real = 0.05, beta::Real = 0.2,\n  minconditionalpower::Real = 0.0, MONOTONECONDITIONALPOWER::Bool = true,\n  npriorpivots::Integer = 50, ngpivots::Integer = 15,\n  label::String = \"\"\n)\n\nThis type represents a set of parameters for finding optimal two-stage designs minimizing expected sample size subject to a constraint on expected (transformed) power.\n\nParameters\n\nParameter Description\nsamplespace a sample space object\np0 upper boundary of the null hypothesis\npmcrv point alternative to power on\nprior prior distribution on true response rate, expected score under the prior conditional on a relevant effect is minimized\na scaling parameter of nonlinear power transfrom (set to infinity for thresholding)\nb scaling parameter of nonlinear power transfrom\ntargetpower target power for thresholding\nk factor for non-responder costs\nalpha maximal type one error rate\nbeta expected transformed power must be larger than 1 - beta\nminconditionalpower minimal conditional power upon continuation to stage two\nMONOTONECONDITIONALPOWER if true, the conditional power must be monotonously increasing, this constraint is only relevant if nmax is set very restrictively\nnpriorpivots number of pivot points for integration of prior\nnpivots number of pivot points piecewise linear approximation of local score\nlabel descriptive string for parameter set\n\n\n\n"
},

{
    "location": "parameters.html#MBESS-1",
    "page": "Design Parameters",
    "title": "MBESS",
    "category": "section",
    "text": "MBESS"
},

{
    "location": "designs.html#",
    "page": "Two-Stage Designs",
    "title": "Two-Stage Designs",
    "category": "page",
    "text": ""
},

{
    "location": "designs.html#BinaryTwoStageDesigns.Design",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.Design",
    "category": "Type",
    "text": "Design{T1<:Integer, T2<:Real, PType<:Parameters}\n\nBinary two-stage design object defining a discrete sampe size function  n(x_1) and c(x_1), where x_1 is the observed number of stage-one events  and the test rejects the superiority null hypothesis whenever x_1 + x_2  c(x_1).\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.Design-Union{Tuple{Any,Any,PType}, Tuple{PType}} where PType<:BinaryTwoStageDesigns.Parameters",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.Design",
    "category": "Method",
    "text": "Design(n, c, params::PType) where {PType<:Parameters}\n\nConstruct a Design object from n and c directly.  Here, n must be convertable to an integer vector and c to a real vector  (both of same length n_1).  Then, n[x1 + 1] is the final sample size after observing x1 responses in stage one and c[x1 + 1] the corresponding critical value. A constructor without explicitly specifying a set of parameters with respect to which thee design is optimal also exists.\n\nParameters\n\nParameter Description\nn sample size function (must be convertable to integer vector)\nc critical value function (must be convertable to real vector)\nparams parameter object\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.Design-Tuple{Any,Any}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.Design",
    "category": "Method",
    "text": "Design(n, c)\n\nConstruct a Design object from n and c directly.  Here, n must be convertable to an integer vector and c to a real vector  (both of same length n_1). \n\nParameters\n\nParameter Description\nn sample size function (must be convertable to integer vector)\nc critical value function (must be convertable to real vector)\n\n\n\n"
},

{
    "location": "designs.html#Base.convert-Tuple{Type{DataFrames.DataFrame},BinaryTwoStageDesigns.Design}",
    "page": "Two-Stage Designs",
    "title": "Base.convert",
    "category": "Method",
    "text": "convert(::Type{DataFrames.DataFrame}, design::Design)\n\nConvert a design to a DataFrame with column names x1, n, and c.\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.interimsamplesize-Tuple{BinaryTwoStageDesigns.Design}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.interimsamplesize",
    "category": "Method",
    "text": "interimsamplesize(design::Design)\n\nReturn interim sample size n_1 of a design.\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.parameters-Tuple{BinaryTwoStageDesigns.Design}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.parameters",
    "category": "Method",
    "text": "parameters(design::Design)\n\nReturn stored parameter object of a design.\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.samplesize-Tuple{BinaryTwoStageDesigns.Design}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.samplesize",
    "category": "Method",
    "text": "samplesize(design::Design)\n\nReturn vector of final sample sizes.\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.samplesize-Union{Tuple{BinaryTwoStageDesigns.Design,T}, Tuple{T}} where T<:Integer",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.samplesize",
    "category": "Method",
    "text": "samplesize(design::Design, x1::T) where {T<:Integer}\n\nReturn final sample size for x1.\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.criticalvalue-Tuple{BinaryTwoStageDesigns.Design}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.criticalvalue",
    "category": "Method",
    "text": "criticalvalue(design::Design)\n\nReturn vector of critical values.\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.criticalvalue-Union{Tuple{BinaryTwoStageDesigns.Design,T}, Tuple{T}} where T<:Integer",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.criticalvalue",
    "category": "Method",
    "text": "criticalvalue(design::Design, x1::T) where {T<:Integer}\n\nReturn critical value for x1.\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.power-Union{Tuple{BinaryTwoStageDesigns.Design,T}, Tuple{T}} where T<:Real",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.power",
    "category": "Method",
    "text": "power(design::Design, p::T) where {T<:Real}\n\nPower of a design for a given response rate p .\n\nParameters\n\nParameter Description\ndesign a Design\np response probability\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.power-Union{Tuple{BinaryTwoStageDesigns.Design,T1,T2}, Tuple{T1}, Tuple{T2}} where T2<:Real where T1<:Integer",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.power",
    "category": "Method",
    "text": "power(design::Design, x1::T1, p::T2) where {T1<:Integer, T2<:Real}\n\nconditional power of a design for a given response rate p and the observed number of stage-one responses x1.\n\nParameters\n\nParameter Description\ndesign a Design\np response probability\nx1 number of stage-one responses\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.expectedpower-Union{Tuple{BinaryTwoStageDesigns.Design,T,Function}, Tuple{T}} where T<:Integer",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.expectedpower",
    "category": "Method",
    "text": "expectedpower(\n  design::Design, x1::T, prior::Function; mcrv::Real = mcrv(parameters(design))\n) where {T<:Integer}\n\nCompute the conditional expected power of a design given x1.\n\nParameters\n\nParameter Description\ndesign a Design\nx1 stage one number of responses\nprior prior function prior(p) for response probability p\nmcrv minimal clinically relevant value for expected power calculation\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.expectedpower-Tuple{BinaryTwoStageDesigns.Design,Function}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.expectedpower",
    "category": "Method",
    "text": "expectedpower(design::Design, prior::Function; mcrv::Real = mcrv(parameters(design)))\n\nCompute the expected power of a given design.\n\nParameters\n\nParameter Description\ndesign a Design\nprior prior function prior(p) for response probability p\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.stoppingforfutility-Union{Tuple{BinaryTwoStageDesigns.Design,T}, Tuple{T}} where T<:Real",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.stoppingforfutility",
    "category": "Method",
    "text": "stoppingforfutility{T<:Real}(design::Design, p::T) where {T<:Real}\n\nCompute probability of stopping early for futility of a given design.\n\nParameters\n\nParameter Description\ndesign a Design\np response probability\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.stoppingforefficacy-Union{Tuple{BinaryTwoStageDesigns.Design,T}, Tuple{T}} where T<:Real",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.stoppingforefficacy",
    "category": "Method",
    "text": "stoppingforefficacy{T<:Real}(design::Design, p::T) where {T<:Real}\n\nCompute probability of stopping early for futility of a given design.\n\nParameters\n\nParameter Description\ndesign a Design\np response probability\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.score-Tuple{BinaryTwoStageDesigns.Design,BinaryTwoStageDesigns.Parameters}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.score",
    "category": "Method",
    "text": "score(design::Design, params::Parameters)::Real\n\nEvaluates the score of a design for a given set of parameters.\n\nParameters\n\nParameter Description\ndesign a Design\nparams a parameters object holding the required information about the score\n\nDetails\n\nA method score(design::Design)::Real is also implemented which uses the parameter object stored within the design object after optimization. Note that this is only possible if the design was created as result of calling getoptimaldesign().\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.test-Union{Tuple{BinaryTwoStageDesigns.Design,T,T}, Tuple{T}} where T<:Integer",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.test",
    "category": "Method",
    "text": "test(design::Design, x1::T, x2::T)::Bool where {T<:Integer}\n\nBinary test decision of design when observing x1 responses in stage one and x2 responses in stage two.\n\nParameters\n\nParameter Description\ndesign a Design\nx1 number of stage-one responses\nx2 number of stage-two responses\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.simulate-Union{Tuple{BinaryTwoStageDesigns.Design,T2,T1}, Tuple{T1}, Tuple{T2}} where T2<:Real where T1<:Integer",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.simulate",
    "category": "Method",
    "text": "simulate(design::Design, p::T2, nsim::T1) where {T1<:Integer, T2<:Real}\n\nSimulate nsim trials of design design with true response probability p.\n\nParameters\n\nParameter Description\ndesign a Design\np true response probability\nnsim number of simulated trials\n\nReturn Value\n\nA DataFrame with columns x1 (stage one responses) n (overall sample size) c (critical value), x2 (stage-two responses), and rejectedH0 (test decision).\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.jeffreysprior-Tuple{BinaryTwoStageDesigns.Design}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.jeffreysprior",
    "category": "Method",
    "text": "jeffreysprior(design::Design)\n\nComputes the Jeffreys prior of any given design.\n\nParameters\n\nParameter Description\ndesign a Design\n\nReturn Value\n\nA function with signature prior{T<:Real}(p::T)::Real where p is the response probability and prior(p) the PDF of the Jeffres prior at p.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> params = SimpleMinimalExpectedSampleSize(ss, .2, .4, .05, .2, .4)\njulia> design, res = getoptimaldesign(params, solver = Gurobi.GurobiSolver())\njulia> f = jeffreysprior(design)\njulia> f(.5)\n\n\n\n"
},

{
    "location": "designs.html#Distributions.pdf-Union{Tuple{BinaryTwoStageDesigns.Design,T1,T1,T2}, Tuple{T1}, Tuple{T2}} where T2<:Real where T1<:Integer",
    "page": "Two-Stage Designs",
    "title": "Distributions.pdf",
    "category": "Method",
    "text": "pdf{T1<:Integer, T2<:Real}(design::Design, x1::T1, x2::T1, p::T2)\n\nProbability density function of (x1 x2) responses in stage one and two,  respectively, under design given response rate p.\n\nParameters\n\nParameter Description\ndesign a Design\nx1 number of stage-one responses\nx2 number of stage-two responses\np response probability\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.save-Tuple{String,BinaryTwoStageDesigns.Design}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.save",
    "category": "Method",
    "text": "save(filename::String, design::Design)\n\nSave design as .jls file.\n\nParameters\n\nParameter Description\nfilename filename for sving design\ndesign a Design\n\n\n\n"
},

{
    "location": "designs.html#Base.DataFmt.writecsv-Tuple{String,BinaryTwoStageDesigns.Design}",
    "page": "Two-Stage Designs",
    "title": "Base.DataFmt.writecsv",
    "category": "Method",
    "text": "writecsv(filename::String, design::Design; label::String = \"\")\n\nSave design as .csv file.\n\nParameters\n\nParameter Description\nfilename filename for sving design\ndesign a Design\nlabel string name of design\n\n\n\n"
},

{
    "location": "designs.html#Binary-Two-Stage-Designs-1",
    "page": "Two-Stage Designs",
    "title": "Binary Two-Stage Designs",
    "category": "section",
    "text": "Design{T1<:Integer, T2<:Real, PType<:Parameters}\n\nDesign(n, c, params::PType) where {PType<:Parameters}\n\nDesign(n, c)\n\nconvert(::Type{DataFrames.DataFrame}, design::Design)\n\ninterimsamplesize(design::Design)\n\nparameters(design::Design)\n\nsamplesize(design::Design)\n\nsamplesize(design::Design, x1::T) where {T<:Integer}\n\ncriticalvalue(design::Design)\n\ncriticalvalue(design::Design, x1::T) where {T<:Integer}\n\npower{T<:Real}(design::Design, p::T)\n\npower(design::Design, x1::T1, p::T2) where {T1<:Integer, T2<:Real}\n\nexpectedpower(design::Design, x1::T, prior::Function; mcrv::Real = mcrv(parameters(design))) where {T<:Integer}\n\nexpectedpower(design::Design, prior::Function; mcrv::Real = mcrv(parameters(design)))\n\nstoppingforfutility{T<:Real}(design::Design, p::T)\n\nstoppingforefficacy{T<:Real}(design::Design, p::T)\n\nscore(design::Design, params::Parameters)\n\ntest(design::Design, x1::T, x2::T) where {T<:Integer}\n\nsimulate(design::Design, p::T2, nsim::T1) where {T1<:Integer, T2<:Real}\n\njeffreysprior(design::Design)\n\npdf(design::Design, x1::T1, x2::T1, p::T2) where {T1<:Integer, T2<:Real}\n\nsave(filename::String, design::Design)\n\nwritecsv(filename::String, design::Design; label::String = \"\") "
},

{
    "location": "optimal_designs.html#",
    "page": "Finding Optimal Designs",
    "title": "Finding Optimal Designs",
    "category": "page",
    "text": ""
},

{
    "location": "optimal_designs.html#BinaryTwoStageDesigns.optimaldesign-Tuple{Integer,BinaryTwoStageDesigns.Parameters,MathProgBase.SolverInterface.AbstractMathProgSolver}",
    "page": "Finding Optimal Designs",
    "title": "BinaryTwoStageDesigns.optimaldesign",
    "category": "Method",
    "text": "optimaldesign(\n  n1::Integer,\n  parameters::Parameters,\n  solver::MathProgBase.AbstractMathProgSolver\n)\n\nFind the optimal two-stage design for given parameters and fixed stage-one sample size.\n\nParameters\n\nParameter Description\nn1 stage-one sample size\nnparameters paramters object defining the optimality criterion\nsolver MathProgBase solver used for optimization\nVERBOSE control verbosity during optimization\n\nReturn Value\n\nDesign object optimized for given parameter set.\n\nExamples\n\njulia> ss     = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> design = optimaldesign(15, params, solver = Gurobi.GurobiSolver())\n\n\n\n"
},

{
    "location": "optimal_designs.html#Optimal-Two-Stage-Designs-1",
    "page": "Finding Optimal Designs",
    "title": "Optimal Two-Stage Designs",
    "category": "section",
    "text": "The technical background in available at arxiv.org.optimaldesign(n1::Integer, parameters::Parameters, solver::MathProgBase.AbstractMathProgSolver)"
},

]}
