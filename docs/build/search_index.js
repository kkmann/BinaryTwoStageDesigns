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
    "text": "Depth = 3\r\nPages = [\r\n  \"sample_space.md\",\r\n  \"parameters.md\",\r\n  \"designs.md\",\r\n  \"optimal_designs.md\"\r\n]"
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
    "category": "type",
    "text": "SampleSpace\n\nDefines restrictions on the sample space of a two-stage design to be used during optimization.\n\n\n\n"
},

{
    "location": "sample_space.html#BinaryTwoStageDesigns.SampleSpace-Union{Tuple{Any,TI}, Tuple{TI}, Tuple{TR}} where TR<:Real where TI<:Integer",
    "page": "Sample Spaces",
    "title": "BinaryTwoStageDesigns.SampleSpace",
    "category": "method",
    "text": "SampleSpace(\n  n1range,\n  nmax::TI;\n  n2min::TI = 1, \n  maxnfact::TR = Inf, \n  nmincont::TI = 0, \n  maxvariables::TI = 500000, \n  GS::Bool = false\n) where {TI<:Integer,TR<:Real}\n\nCreates a SampleSpace-object which defines the search space for finding  optimal two-stage designs using optimaldesign. The primary use of dedicated sample sapce objects is the incoporation of  specific operational constraints on optimal designs (e.g. maximal sample size).\n\nParameters\n\nParameter Description\nn1range possible stage-one sample sizes, can be anything which is convertable to an integer vector\nnmax maximal overall sample size (stage one and two combined)\n[n2min] minimal stage-two sample size\n[maxnfact] maximal overall sample size must be smaller than maxnfact * n1\n[nmincont] minimal sample size upon rejection of the null hypothesis\n[maxvariables] (approximate) maximal number of binary variables used when finding optimal design; if the sample space is too big this can be used to find approximate solutions to the optimization problem\n[GS] flag indicating wheather the design should be group-sequential (constant stage two sample size upon continuation)\n\n\n\n"
},

{
    "location": "sample_space.html#BinaryTwoStageDesigns.interimsamplesizerange-Tuple{BinaryTwoStageDesigns.SampleSpace}",
    "page": "Sample Spaces",
    "title": "BinaryTwoStageDesigns.interimsamplesizerange",
    "category": "method",
    "text": "interimsamplesizerange(ss::SampleSpace)\n\nReturn integer vector of valid n_1 values.\n\n\n\n"
},

{
    "location": "sample_space.html#BinaryTwoStageDesigns.maxsamplesize-Tuple{BinaryTwoStageDesigns.SampleSpace}",
    "page": "Sample Spaces",
    "title": "BinaryTwoStageDesigns.maxsamplesize",
    "category": "method",
    "text": "maxsamplesize(ss::SampleSpace)\n\nReturn maximal overal sample size as integer.\n\n\n\n"
},

{
    "location": "sample_space.html#BinaryTwoStageDesigns.maxsamplesize-Union{Tuple{BinaryTwoStageDesigns.SampleSpace{TI,TR},TI2}, Tuple{TI2}, Tuple{TI}, Tuple{TR}} where TI2<:Integer where TR<:Real where TI<:Integer",
    "page": "Sample Spaces",
    "title": "BinaryTwoStageDesigns.maxsamplesize",
    "category": "method",
    "text": "maxsamplesize(\n  ss::SampleSpace{TI,TR}, n1::TI2\n) where {TI<:Integer,TR<:Real,TI2<:Integer}\n\nReturn maximal sample size given n1 as integer.\n\n\n\n"
},

{
    "location": "sample_space.html#BinaryTwoStageDesigns.isgroupsequential-Tuple{BinaryTwoStageDesigns.SampleSpace}",
    "page": "Sample Spaces",
    "title": "BinaryTwoStageDesigns.isgroupsequential",
    "category": "method",
    "text": "isgroupsequential(ss::SampleSpace)\n\nIndicate whether the sample space is restricted to group-sequential designs.\n\n\n\n"
},

{
    "location": "sample_space.html#The-Sample-Space-1",
    "page": "Sample Spaces",
    "title": "The Sample Space",
    "category": "section",
    "text": "Sample space objects can be used to encode the feasible space of binary  two-stage designs for optimization.SampleSpace\r\n\r\n SampleSpace(\r\n  n1range,\r\n  nmax::TI;\r\n  n2min::TI = 1, \r\n  maxnfact::TR = Inf, \r\n  nmincont::TI = 0, \r\n  maxvariables::TI = 500000, \r\n  GS::Bool = false\r\n) where {TI<:Integer,TR<:Real}\r\n\r\ninterimsamplesizerange(ss::SampleSpace)\r\n\r\nmaxsamplesize(ss::SampleSpace)\r\n\r\nmaxsamplesize(ss::SampleSpace{TI,TR}, n1::TI2) where {TI<:Integer,TR<:Real,TI2<:Integer}\r\n\r\nisgroupsequential(ss::SampleSpace)"
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
    "category": "type",
    "text": "Parameters\n\nAbstract type representing a generic set of paramters for finding optimal two-stage designs. \n\n\n\n"
},

{
    "location": "parameters.html#BinaryTwoStageDesigns.label-Tuple{BinaryTwoStageDesigns.Parameters}",
    "page": "Design Parameters",
    "title": "BinaryTwoStageDesigns.label",
    "category": "method",
    "text": "label(par::Parameters)\n\nReturn the label of a Parameters-object.\n\n\n\n"
},

{
    "location": "parameters.html#BinaryTwoStageDesigns.null-Tuple{BinaryTwoStageDesigns.Parameters}",
    "page": "Design Parameters",
    "title": "BinaryTwoStageDesigns.null",
    "category": "method",
    "text": "null(par::Parameters)\n\nReturn the response rate under the null hypothesis of a Parameters-object.\n\n\n\n"
},

{
    "location": "parameters.html#BinaryTwoStageDesigns.mtoer-Tuple{BinaryTwoStageDesigns.Parameters}",
    "page": "Design Parameters",
    "title": "BinaryTwoStageDesigns.mtoer",
    "category": "method",
    "text": "mtoer(par::Parameters)\n\nReturn the maximal type one error rate of a Parameters-object.\n\n\n\n"
},

{
    "location": "parameters.html#BinaryTwoStageDesigns.mcrv-Tuple{BinaryTwoStageDesigns.Parameters}",
    "page": "Design Parameters",
    "title": "BinaryTwoStageDesigns.mcrv",
    "category": "method",
    "text": "mcrv(par::Parameters)\n\nReturn the minimal clinically relevant response rate of a Parameters-object.\n\n\n\n"
},

{
    "location": "parameters.html#BinaryTwoStageDesigns.samplespace-Tuple{BinaryTwoStageDesigns.Parameters}",
    "page": "Design Parameters",
    "title": "BinaryTwoStageDesigns.samplespace",
    "category": "method",
    "text": "samplespace(par::Parameters)\n\nReturn the SampleSpace-object of a Parameters-object.\n\n\n\n"
},

{
    "location": "parameters.html#Overview-1",
    "page": "Design Parameters",
    "title": "Overview",
    "category": "section",
    "text": "Objects of type Parameters are used to define objective functions and their  respective parameters for finding optimal two-stage designs. Every Parameters object also contains a SampleSpace which restricts  the feasible region of the sample space and implements a score method  which is to be optimized.Parameters\r\n\r\nlabel(par::Parameters)\r\n\r\nnull(par::Parameters)\r\n\r\nmtoer(par::Parameters)\r\n\r\nmcrv(par::Parameters)\r\n\r\nsamplespace(par::Parameters)"
},

{
    "location": "parameters.html#BinaryTwoStageDesigns.MESS",
    "page": "Design Parameters",
    "title": "BinaryTwoStageDesigns.MESS",
    "category": "type",
    "text": "MESS{T_samplespace<:SampleSpace,TR<:Real} <: PointAlternative\n\nMESS{T_samplespace<:SampleSpace}(\n  samplespace::T_samplespace,\n  p0, p1,\n  alpha, beta,\n  pess;\n  minstoppingforfutility::Real   = 0.0,\n  minconditionalpower::Real      = 0.0,\n  MONOTONECONDITIONALPOWER::Bool = true\n)\n\nThis type represents a set of parameters for finding optimal two-stage designs minimizing the expected sample size on a point in the parameter space subject to type one and two error rate constraints (Minimal Expected Sample Size).\n\n> [1] Simon R. Optimal two-stage designs for phase II clinical trials. `Controlled Clinical Trials` 1989; 10, 1-10.\n\n> [2] Kunzmann K and Kieser M. Optimal adaptive two-stage designs for single-arm trials with binary endpoint. `arxive.org` 2016; arXiv:1605.00249.\n\nDetails\n\nPlease note that due to technical reasons the power constraint is actually imlemented in a relaxed way by adding lambda*(1-beta - power)_+ to the objective function. In practice, lambda should be large enough to ensure a power of 1-beta, but depending on the setting, a larger lambda might be required.\n\nParameters\n\nParameter Description\nsamplespace a sample space object\np0 upper boundary of the null hypothesis\np1 point alternative to power on\nalpha maximal tolerable type one error rate\nbeta maximal tolerable type two error rate on p1\npess response rate under which to minimize expected sample size\nminstoppingforfutility minimal probability for stopping for futility under p0\nminconditionalpower minimal conditional power upon continuation to stage two\nMONOTONECONDITIONALPOWER if true, the conditional power must be monotonously increasing, this constraint is only relevant if nmax is set very restrictively\n\n\n\n"
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
    "category": "type",
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
    "category": "type",
    "text": "EB{TI<:Integer,TR<:Real} <: VagueAlternative\n\nEB( # default values\n  samplespace::SampleSpace,\n  p0::Real, pmcrv::Real, prior,\n  gamma::Real, lambda::Real;\n  a::Real = 1, b::Real = 1, targetpower::Real = .8, k::Real = 1,\n  alpha::Real = 0.05,\n  minconditionalpower::Real = 0.0, MONOTONECONDITIONALPOWER::Bool = true,\n  npriorpivots::Integer = 50, ngpivots::Integer = 15,\n  label::String = \"\"\n)\n\nThis type represents a set of parameters for finding optimal two-stage designs maximizing the expected benefit (EB).\n\nParameters\n\nParameter Description\nsamplespace a sample space object\np0 upper boundary of the null hypothesis\npmcrv point alternative to power on\nprior prior distribution on true response rate, expected score under the prior conditional on a relevant effect is minimized\ngamma per patient costs\nlambda benefit of correctly rejecting null if relevant effect exists\na scaling parameter of nonlinear power transfrom (set to infinity for thresholding)\nb scaling parameter of nonlinear power transfrom\ntargetpower target power for thresholding\nk factor for non-responder costs\nalpha maximal type one error rate\nminconditionalpower minimal conditional power upon continuation to stage two\nMONOTONECONDITIONALPOWER if true, the conditional power must be monotonously increasing, this constraint is only relevant if nmax is set very restrictively\nnpriorpivots number of pivot points for integration of prior\nnpivots number of pivot points piecewise linear approximation of local score\nlabel descriptive string for parameter set\n\n\n\n"
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
    "category": "type",
    "text": "MBESS{TI<:Integer,TR<:Real} <: VagueAlternative\n\nMBESS(\n  samplespace::SampleSpace,\n  p0::Real, pmcrv::Real, prior;\n  a::Real = 1, b::Real = 1, targetpower::Real = .8, k::Real = 1,\n  alpha::Real = 0.05, beta::Real = 0.2,\n  minconditionalpower::Real = 0.0, MONOTONECONDITIONALPOWER::Bool = true,\n  npriorpivots::Integer = 50, ngpivots::Integer = 15,\n  label::String = \"\"\n)\n\nThis type represents a set of parameters for finding optimal two-stage designs minimizing expected sample size subject to a constraint on expected (transformed) power (Minimal Bayesian Expected Sample Size).\n\nParameters\n\nParameter Description\nsamplespace a sample space object\np0 upper boundary of the null hypothesis\npmcrv point alternative to power on\nprior prior distribution on true response rate, expected score under the prior conditional on a relevant effect is minimized\na scaling parameter of nonlinear power transfrom (set to infinity for thresholding)\nb scaling parameter of nonlinear power transfrom\ntargetpower target power for thresholding\nk factor for non-responder costs\nalpha maximal type one error rate\nbeta expected transformed power must be larger than 1 - beta\nminconditionalpower minimal conditional power upon continuation to stage two\nMONOTONECONDITIONALPOWER if true, the conditional power must be monotonously increasing, this constraint is only relevant if nmax is set very restrictively\nnpriorpivots number of pivot points for integration of prior\nnpivots number of pivot points piecewise linear approximation of local score\nlabel descriptive string for parameter set\n\n\n\n"
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
    "category": "type",
    "text": "Design{T1<:Integer, T2<:Real, PType<:Parameters}\n\nBinary two-stage design object defining a discrete sampe size function  n(x_1) and c(x_1), where x_1 is the observed number of stage-one events  and the test rejects the superiority null hypothesis whenever x_1 + x_2  c(x_1).\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.Design-Union{Tuple{Any,Any,PType}, Tuple{PType}} where PType<:BinaryTwoStageDesigns.Parameters",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.Design",
    "category": "method",
    "text": "Design(n, c, params::PType) where {PType<:Parameters}\n\nConstruct a Design object from n and c directly.  Here, n must be convertable to an integer vector and c to a real vector  (both of same length n_1).  Then, n[x1 + 1] is the final sample size after observing x1 responses in stage one and c[x1 + 1] the corresponding critical value. A constructor without explicitly specifying a set of parameters with respect to which thee design is optimal also exists.\n\nParameters\n\nParameter Description\nn sample size function (must be convertable to integer vector)\nc critical value function (must be convertable to real vector)\nparams parameter object\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.Design-Tuple{Any,Any}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.Design",
    "category": "method",
    "text": "Design(n, c)\n\nConstruct a Design object from n and c directly.  Here, n must be convertable to an integer vector and c to a real vector  (both of same length n_1). \n\nParameters\n\nParameter Description\nn sample size function (must be convertable to integer vector)\nc critical value function (must be convertable to real vector)\n\n\n\n"
},

{
    "location": "designs.html#Base.convert-Tuple{Type{DataFrames.DataFrame},BinaryTwoStageDesigns.Design}",
    "page": "Two-Stage Designs",
    "title": "Base.convert",
    "category": "method",
    "text": "convert(::Type{DataFrames.DataFrame}, design::Design)\n\nConvert a design to a DataFrame with column names x1, n, and c.\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.interimsamplesize-Tuple{BinaryTwoStageDesigns.Design}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.interimsamplesize",
    "category": "method",
    "text": "interimsamplesize(design::Design)\n\nReturn interim sample size n_1 of a design.\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.parameters-Tuple{BinaryTwoStageDesigns.Design}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.parameters",
    "category": "method",
    "text": "parameters(design::Design)\n\nReturn stored parameter object of a design.\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.samplesize-Tuple{BinaryTwoStageDesigns.Design}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.samplesize",
    "category": "method",
    "text": "samplesize(design::Design)\n\nReturn vector of final sample sizes.\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.samplesize-Union{Tuple{BinaryTwoStageDesigns.Design,T}, Tuple{T}} where T<:Integer",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.samplesize",
    "category": "method",
    "text": "samplesize(design::Design, x1::T) where {T<:Integer}\n\nReturn final sample size for x1.\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.criticalvalue-Tuple{BinaryTwoStageDesigns.Design}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.criticalvalue",
    "category": "method",
    "text": "criticalvalue(design::Design)\n\nReturn vector of critical values.\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.criticalvalue-Union{Tuple{BinaryTwoStageDesigns.Design,T}, Tuple{T}} where T<:Integer",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.criticalvalue",
    "category": "method",
    "text": "criticalvalue(design::Design, x1::T) where {T<:Integer}\n\nReturn critical value for x1.\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.power-Union{Tuple{BinaryTwoStageDesigns.Design,T}, Tuple{T}} where T<:Real",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.power",
    "category": "method",
    "text": "power(design::Design, p::T) where {T<:Real}\n\nPower of a design for a given response rate p .\n\nParameters\n\nParameter Description\ndesign a Design\np response probability\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.power-Union{Tuple{BinaryTwoStageDesigns.Design,T1,T2}, Tuple{T1}, Tuple{T2}} where T2<:Real where T1<:Integer",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.power",
    "category": "method",
    "text": "power(design::Design, x1::T1, p::T2) where {T1<:Integer, T2<:Real}\n\nconditional power of a design for a given response rate p and the observed number of stage-one responses x1.\n\nParameters\n\nParameter Description\ndesign a Design\np response probability\nx1 number of stage-one responses\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.expectedpower-Union{Tuple{BinaryTwoStageDesigns.Design,T,Function}, Tuple{T}} where T<:Integer",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.expectedpower",
    "category": "method",
    "text": "expectedpower(\n  design::Design, x1::T, prior::Function; mcrv::Real = mcrv(parameters(design))\n) where {T<:Integer}\n\nCompute the conditional expected power of a design given x1.\n\nParameters\n\nParameter Description\ndesign a Design\nx1 stage one number of responses\nprior prior function prior(p) for response probability p\nmcrv minimal clinically relevant value for expected power calculation\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.expectedpower-Tuple{BinaryTwoStageDesigns.Design,Function}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.expectedpower",
    "category": "method",
    "text": "expectedpower(design::Design, prior::Function; mcrv::Real = mcrv(parameters(design)))\n\nCompute the expected power of a given design.\n\nParameters\n\nParameter Description\ndesign a Design\nprior prior function prior(p) for response probability p\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.stoppingforfutility-Union{Tuple{BinaryTwoStageDesigns.Design,T}, Tuple{T}} where T<:Real",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.stoppingforfutility",
    "category": "method",
    "text": "stoppingforfutility{T<:Real}(design::Design, p::T) where {T<:Real}\n\nCompute probability of stopping early for futility of a given design.\n\nParameters\n\nParameter Description\ndesign a Design\np response probability\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.stoppingforefficacy-Union{Tuple{BinaryTwoStageDesigns.Design,T}, Tuple{T}} where T<:Real",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.stoppingforefficacy",
    "category": "method",
    "text": "stoppingforefficacy{T<:Real}(design::Design, p::T) where {T<:Real}\n\nCompute probability of stopping early for futility of a given design.\n\nParameters\n\nParameter Description\ndesign a Design\np response probability\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.score-Tuple{BinaryTwoStageDesigns.Design,BinaryTwoStageDesigns.Parameters}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.score",
    "category": "method",
    "text": "score(design::Design, params::Parameters)::Real\n\nEvaluates the score of a design for a given set of parameters.\n\nParameters\n\nParameter Description\ndesign a Design\nparams a parameters object holding the required information about the score\n\nDetails\n\nA method score(design::Design)::Real is also implemented which uses the parameter object stored within the design object after optimization. Note that this is only possible if the design was created as result of calling getoptimaldesign().\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.test-Union{Tuple{BinaryTwoStageDesigns.Design,T,T}, Tuple{T}} where T<:Integer",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.test",
    "category": "method",
    "text": "test(design::Design, x1::T, x2::T)::Bool where {T<:Integer}\n\nBinary test decision of design when observing x1 responses in stage one and x2 responses in stage two.\n\nParameters\n\nParameter Description\ndesign a Design\nx1 number of stage-one responses\nx2 number of stage-two responses\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.simulate-Union{Tuple{BinaryTwoStageDesigns.Design,T2,T1}, Tuple{T1}, Tuple{T2}} where T2<:Real where T1<:Integer",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.simulate",
    "category": "method",
    "text": "simulate(design::Design, p::T2, nsim::T1) where {T1<:Integer, T2<:Real}\n\nSimulate nsim trials of design design with true response probability p.\n\nParameters\n\nParameter Description\ndesign a Design\np true response probability\nnsim number of simulated trials\n\nReturn Value\n\nA DataFrame with columns x1 (stage one responses) n (overall sample size) c (critical value), x2 (stage-two responses), and rejectedH0 (test decision).\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.jeffreysprior-Tuple{BinaryTwoStageDesigns.Design}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.jeffreysprior",
    "category": "method",
    "text": "jeffreysprior(design::Design)\n\nComputes the Jeffreys prior of any given design.\n\nParameters\n\nParameter Description\ndesign a Design\n\nReturn Value\n\nA function with signature prior{T<:Real}(p::T)::Real where p is the response probability and prior(p) the PDF of the Jeffres prior at p.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> params = SimpleMinimalExpectedSampleSize(ss, .2, .4, .05, .2, .4)\njulia> design, res = getoptimaldesign(params, solver = Gurobi.GurobiSolver())\njulia> f = jeffreysprior(design)\njulia> f(.5)\n\n\n\n"
},

{
    "location": "designs.html#Distributions.pdf-Union{Tuple{BinaryTwoStageDesigns.Design,T1,T1,T2}, Tuple{T1}, Tuple{T2}} where T2<:Real where T1<:Integer",
    "page": "Two-Stage Designs",
    "title": "Distributions.pdf",
    "category": "method",
    "text": "pdf{T1<:Integer, T2<:Real}(design::Design, x1::T1, x2::T1, p::T2)\n\nProbability density function of (x1 x2) responses in stage one and two,  respectively, under design given response rate p.\n\nParameters\n\nParameter Description\ndesign a Design\nx1 number of stage-one responses\nx2 number of stage-two responses\np response probability\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.save-Tuple{String,BinaryTwoStageDesigns.Design}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.save",
    "category": "method",
    "text": "save(filename::String, design::Design)\n\nSave design as .jls file.\n\nParameters\n\nParameter Description\nfilename filename for sving design\ndesign a Design\n\n\n\n"
},

{
    "location": "designs.html#Base.DataFmt.writecsv-Tuple{String,BinaryTwoStageDesigns.Design}",
    "page": "Two-Stage Designs",
    "title": "Base.DataFmt.writecsv",
    "category": "method",
    "text": "writecsv(filename::String, design::Design; label::String = \"\")\n\nSave design as .csv file.\n\nParameters\n\nParameter Description\nfilename filename for sving design\ndesign a Design\nlabel string name of design\n\n\n\n"
},

{
    "location": "designs.html#Binary-Two-Stage-Designs-1",
    "page": "Two-Stage Designs",
    "title": "Binary Two-Stage Designs",
    "category": "section",
    "text": "Design{T1<:Integer, T2<:Real, PType<:Parameters}\r\n\r\nDesign(n, c, params::PType) where {PType<:Parameters}\r\n\r\nDesign(n, c)\r\n\r\nconvert(::Type{DataFrames.DataFrame}, design::Design)\r\n\r\ninterimsamplesize(design::Design)\r\n\r\nparameters(design::Design)\r\n\r\nsamplesize(design::Design)\r\n\r\nsamplesize(design::Design, x1::T) where {T<:Integer}\r\n\r\ncriticalvalue(design::Design)\r\n\r\ncriticalvalue(design::Design, x1::T) where {T<:Integer}\r\n\r\npower{T<:Real}(design::Design, p::T)\r\n\r\npower(design::Design, x1::T1, p::T2) where {T1<:Integer, T2<:Real}\r\n\r\nexpectedpower(design::Design, x1::T, prior::Function; mcrv::Real = mcrv(parameters(design))) where {T<:Integer}\r\n\r\nexpectedpower(design::Design, prior::Function; mcrv::Real = mcrv(parameters(design)))\r\n\r\nstoppingforfutility{T<:Real}(design::Design, p::T)\r\n\r\nstoppingforefficacy{T<:Real}(design::Design, p::T)\r\n\r\nscore(design::Design, params::Parameters)\r\n\r\ntest(design::Design, x1::T, x2::T) where {T<:Integer}\r\n\r\nsimulate(design::Design, p::T2, nsim::T1) where {T1<:Integer, T2<:Real}\r\n\r\njeffreysprior(design::Design)\r\n\r\npdf(design::Design, x1::T1, x2::T1, p::T2) where {T1<:Integer, T2<:Real}\r\n\r\nsave(filename::String, design::Design)\r\n\r\nwritecsv(filename::String, design::Design; label::String = \"\") "
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
    "category": "method",
    "text": "optimaldesign(\n  n1::Integer,\n  parameters::Parameters,\n  solver::MathProgBase.AbstractMathProgSolver\n)\n\nFind the optimal two-stage design for given parameters and fixed stage-one sample size.\n\nParameters\n\nParameter Description\nn1 stage-one sample size\nnparameters paramters object defining the optimality criterion\nsolver MathProgBase solver used for optimization\nVERBOSE control verbosity during optimization\n\nReturn Value\n\nDesign object optimized for given parameter set.\n\nExamples\n\njulia> ss     = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> design = optimaldesign(15, params, solver = Gurobi.GurobiSolver())\n\n\n\n"
},

{
    "location": "optimal_designs.html#BinaryTwoStageDesigns.optimaldesign-Tuple{BinaryTwoStageDesigns.Parameters,MathProgBase.SolverInterface.AbstractMathProgSolver}",
    "page": "Finding Optimal Designs",
    "title": "BinaryTwoStageDesigns.optimaldesign",
    "category": "method",
    "text": "getoptimaldesign{TS<:MathProgBase.AbstractMathProgSolver}(\n    parameters::Parameters,\n    solver::TS;\n    VERBOSE::Integer = 1\n)\n\nFind the optimal two-stage design for given parameters (optimizes over n1 as well).\n\nParameters\n\nParameter Description\nnparameters paramters object defining the optimality criterion\nsolver MathProgBase solver used for optimization\nVERBOSE control verbosity during optimization\n\nReturn Value\n\nTuple, first element is the Design object optimized for given parameter set and the second element is a dictionary with n1, scores and respective optimal designs.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design, res = getoptimaldesign(params, solver = Gurobi.GurobiSolver())\n\n\n\n"
},

{
    "location": "optimal_designs.html#Optimal-Two-Stage-Designs-1",
    "page": "Finding Optimal Designs",
    "title": "Optimal Two-Stage Designs",
    "category": "section",
    "text": "The technical background in available at arxiv.org.optimaldesign(n1::Integer, parameters::Parameters, solver::MathProgBase.AbstractMathProgSolver)\r\n\r\noptimaldesign(\r\n  parameters::Parameters,\r\n  solver::MathProgBase.AbstractMathProgSolver;\r\n  VERBOSE::Integer = 1,\r\n  EARLYTERMINATION::Bool = false\r\n)"
},

{
    "location": "inference.html#",
    "page": "Inference",
    "title": "Inference",
    "category": "page",
    "text": ""
},

{
    "location": "inference.html#Inference-1",
    "page": "Inference",
    "title": "Inference",
    "category": "section",
    "text": "Pages = [\"inference.md\"]\r\nDepth = 3"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.Estimator",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.Estimator",
    "category": "type",
    "text": "Estimator\n\nAbstract base type for all estimators.\n\n\n\n"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.estimate-Union{Tuple{BinaryTwoStageDesigns.Estimator,T,T}, Tuple{T}} where T<:Integer",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.estimate",
    "category": "method",
    "text": "estimate{T<:Integer}(estimator::Estimator, x1::T, x2::T)\n\nEstimate the response rate from observed x1 and x2.\n\nParameters\n\nParameter Description\nestimator any BinaryTwoStageDeisgnEstimator object\nx1 stage-one responses\nx2 stage-two responses\n\nReturn Value\n\nReal, estimated response rate.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())\njulia> est = MLE(design)\njulia> estimate(est, 0, 0)\n\n\n\n"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.bias-Union{Tuple{BinaryTwoStageDesigns.Estimator,T}, Tuple{T}} where T<:Real",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.bias",
    "category": "method",
    "text": "bias{T<:Real}(estimator::Estimator, p::T)\n\nBias of estimator given response rate p.\n\nParameters\n\nParameter Description\nestimator any BinaryTwoStageDeisgnEstimator object\np0 upper boundary of null hypothesis\n\nReturn Value\n\nReal, bias given p.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())\njulia> est = MLE(design)\njulia> bias(est, .3)\n\n\n\n"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.rmse-Union{Tuple{BinaryTwoStageDesigns.Estimator,T}, Tuple{T}} where T<:Real",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.rmse",
    "category": "method",
    "text": "rmse{T<:Real}(estimator::Estimator, p::T)\n\nRoot mean squared error of estimator given response rate p.\n\nParameters\n\nParameter Description\nestimator any BinaryTwoStageDeisgnEstimator object\np0 upper boundary of null hypothesis\n\nReturn Value\n\nReal, RMSE given p.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())\njulia> est = MLE(design)\njulia> rmse(est, .3)\n\n\n\n"
},

{
    "location": "inference.html#Point-estimation-1",
    "page": "Inference",
    "title": "Point estimation",
    "category": "section",
    "text": "Estimator\r\n\r\nestimate{T<:Integer}(estimator::Estimator, x1::T, x2::T)\r\n\r\nbias{T<:Real}(estimator::Estimator, p::T)\r\n\r\nrmse{T<:Real}(estimator::Estimator, p::T)"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.MLEstimator-Tuple{BinaryTwoStageDesigns.Design}",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.MLEstimator",
    "category": "method",
    "text": "MLEstimator(design::TD) where {TD<:Design}\n\nCreate maximum likelihood estimator for response rate under given design.\n\n\n\n"
},

{
    "location": "inference.html#Maximum-likelihood-estimator-1",
    "page": "Inference",
    "title": "Maximum likelihood estimator",
    "category": "section",
    "text": "MLEstimator(design::Design)"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.RBEstimator-Tuple{BinaryTwoStageDesigns.Design}",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.RBEstimator",
    "category": "method",
    "text": "RBEstimator(design::TD) where {TD<:Design}\n\nCreate unbiased estimator for response rate p see also:\n\nKunzmann K, Kieser M. Point estimation and p‐values in phase II adaptive two‐stage designs with a binary endpoint. Statistics in medicine. 2017 Mar 15;36(6):971-84.\n\n\n\n"
},

{
    "location": "inference.html#Unbiased-estimator-1",
    "page": "Inference",
    "title": "Unbiased estimator",
    "category": "section",
    "text": "RBEstimator(design::Design)"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.OCEstimator-Union{Tuple{BinaryTwoStageDesigns.Design,TS}, Tuple{TS}} where TS<:MathProgBase.SolverInterface.AbstractMathProgSolver",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.OCEstimator",
    "category": "method",
    "text": "OCEstimator{TS<:MathProgBase.AbstractMathProgSolver}(\n  design::Design,\n  solver::TS;\n  prior::Function = jeffreysprior(design),\n  k = 100\n)\n\nCreate compatible estimator minimizing expected MSE for response rate p see also:\n\nKunzmann K, Kieser M. Point estimation and p‐values in phase II adaptive two‐stage designs with a binary endpoint. Statistics in medicine. 2017 Mar 15;36(6):971-84.\n\nParameters\n\nParameter Description\ndesign Design\nsolver MathProgBase solver, must support quadratic expressions\nprior weight function for MSE values at different p, must be of form f(p::Real)::Real\nk number of equally spaced grid-points for evaluation of MSE and prior\n\n\n\n"
},

{
    "location": "inference.html#Optimal-compatible-estimator-1",
    "page": "Inference",
    "title": "Optimal compatible estimator",
    "category": "section",
    "text": "OCEstimator{TS<:MathProgBase.AbstractMathProgSolver}(\r\n    design::Design,\r\n    solver::TS;\r\n    prior::Function = jeffreysprior(design),\r\n    k = 100\r\n)"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.pvalue-Union{Tuple{BinaryTwoStageDesigns.Estimator,T1,T1,T2}, Tuple{T1}, Tuple{T2}} where T2<:Real where T1<:Integer",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.pvalue",
    "category": "method",
    "text": "p{T1<:Integer, T2<:Real}(estimator::Estimator, x1::T1, x2::T1, p0::T2)\n\nCompute the p value after observing (x1, x2) for null hypothesis H0: p <= p0 with respect to ordering induced by estimator.\n\nParameters\n\nParameter Description\nestimator any BinaryTwoStageDeisgnEstimator object\nx1 stage-one responses\nx2 stage-two responses\np0 upper boundary of null hypothesis\n\nReturn Value\n\nReal, p value.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())\njulia> est = MLE(design)\njulia> p(est, 0, 0, .2)\n\n\n\n"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.incompatibleoutcomes-Tuple{BinaryTwoStageDesigns.Estimator}",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.incompatibleoutcomes",
    "category": "method",
    "text": "incompatibleoutcomes(estimator::Estimator)\n\nFind outcomes where the induced p value implies different decision than the underlying design\n\nParameters\n\nParameter Description\nestimator any BinaryTwoStageDeisgnEstimator object\n\nReturn Value\n\nArray with respective outcomes\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())\njulia> est = MLE(design)\njulia> incompatibleoutcomes(est)\n\n\n\n"
},

{
    "location": "inference.html#P-values-1",
    "page": "Inference",
    "title": "P values",
    "category": "section",
    "text": "pvalue(estimator::Estimator, x1::T1, x2::T1, p0::T2) where {T1<:Integer, T2<:Real}\r\n\r\nincompatibleoutcomes(estimator::Estimator)"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.ConfidenceInterval",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.ConfidenceInterval",
    "category": "type",
    "text": "ConfidenceInterval\n\nAbstract base type for all confidence interval types.\n\n\n\n"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.limits-Union{Tuple{BinaryTwoStageDesigns.ConfidenceInterval,T,T}, Tuple{T}} where T<:Integer",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.limits",
    "category": "method",
    "text": "limits{T<:Integer}(ci::ConfidenceInterval, x1::T, x2::T)\n\nReturn the confidence interval limits for observed x1 and x2.\n\nParameters\n\nParameter Description\nci any ConfidenceInterval object\nx1 stage-one responses\nx2 stage-two responses\n\nReturn Value\n\nTwo element Real vector of limits [lower, upper].\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())\njulia> est = MLE(design)\njulia> ci = ClopperPearsonConfidenceInterval(est, confidence = .9)\njulia> limits(ci, 0, 0)\n\n\n\n"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.coverage-Union{Tuple{BinaryTwoStageDesigns.ConfidenceInterval,T}, Tuple{T}} where T<:Real",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.coverage",
    "category": "method",
    "text": "coverage(ci::ConfidenceInterval, p::T; orientation::String = \"overall\") where {T<:Real}\n\nReturn coverage of given confidence interval and response rate p.\n\nParameters\n\nParameter Description\nci any ConfidenceInterval object\np response rate\norientation string indicating the coverage type - \"overall\", \"lower\", or \"upper\"\n\nReturn Value\n\nCoverage probability given p.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())\njulia> est = MLE(design)\njulia> ci = ClopperPearsonConfidenceInterval(est, confidence = .9)\njulia> coverage(ci, .5)\n\n\n\n"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.meanwidth-Union{Tuple{BinaryTwoStageDesigns.ConfidenceInterval,T}, Tuple{T}} where T<:Real",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.meanwidth",
    "category": "method",
    "text": "meanwidth(ci::ConfidenceInterval, p::T) where {T<:Real}\n\nMean width of given confidence interval and response rate p.\n\nParameters\n\nParameter Description\nci any ConfidenceInterval object\np response rate\n\nReturn Value\n\nMean width given p.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())\njulia> est = MLE(design)\njulia> ci = ClopperPearsonConfidenceInterval(est, confidence = .9)\njulia> meanwidth(ci, .5)\n\n\n\n"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.meaninterval-Union{Tuple{BinaryTwoStageDesigns.ConfidenceInterval,T}, Tuple{T}} where T<:Real",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.meaninterval",
    "category": "method",
    "text": "meaninterval(ci::ConfidenceInterval, p::T) where {T<:Real}\n\nMean interval (average limits) of given confidence interval and response rate p.\n\nParameters\n\nParameter Description\nci any ConfidenceInterval object\np response rate\n\nReturn Value\n\nMean interval given p.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())\njulia> est = MLE(design)\njulia> ci = ClopperPearsonConfidenceInterval(est, confidence = .9)\njulia> meaninterval(ci, .5)\n\n\n\n"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.findinconsistencies-Union{Tuple{BinaryTwoStageDesigns.ConfidenceInterval,T}, Tuple{T}} where T<:Real",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.findinconsistencies",
    "category": "method",
    "text": "findinconsistencies(ci::ConfidenceInterval, p0::T) where {T<:Real}\n\nReturn outcomes where the confidence interval contradicts the designs test decision.\n\nParameters\n\nParameter Description\nci any ConfidenceInterval object\np response rate\n\nReturn Value\n\nMean interval given p.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())\njulia> est = MLE(design)\njulia> ci = ClopperPearsonConfidenceInterval(est, confidence = .9)\njulia> findinconsistencies(ci)\n\n\n\n"
},

{
    "location": "inference.html#Confidence-intervals-1",
    "page": "Inference",
    "title": "Confidence intervals",
    "category": "section",
    "text": "ConfidenceInterval\r\n\r\nlimits(ci::ConfidenceInterval, x1::T, x2::T) where {T<:Integer}\r\n\r\ncoverage(ci::ConfidenceInterval, p::T; orientation::String = \"overall\") where {T<:Real}\r\n\r\nmeanwidth(ci::ConfidenceInterval, p::T) where {T<:Real}\r\n\r\nmeaninterval(ci::ConfidenceInterval, p::T) where {T<:Real}\r\n\r\nfindinconsistencies(ci::ConfidenceInterval, p0::T) where {T<:Real}"
},

{
    "location": "inference.html#Naive-Clopper-Pearson-confidence-interval-1",
    "page": "Inference",
    "title": "Naive Clopper-Pearson confidence interval",
    "category": "section",
    "text": "CPInterval{T<:Real}(\r\n  design::Design;\r\n  confidence::T = .9\r\n)"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.ECPInterval-Union{Tuple{BinaryTwoStageDesigns.Estimator}, Tuple{T}} where T<:Real",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.ECPInterval",
    "category": "method",
    "text": "ECPInterval{T<:Real}(\n  estimator::Estimator;\n  confidence::T = .9,\n  k::Integer = 1001\n)\n\nExact Clopper-Pearson type confidence interval based on ordering induced by estimator.\n\nParameters\n\nParameter Description\nestimator estimator object defining the sample space ordering\nconfidence confidence level of the interval\nk number of equally spaced grid-points for invertign the test\n\n\n\n"
},

{
    "location": "inference.html#Clopper-Pearson-confidence-interval-1",
    "page": "Inference",
    "title": "Clopper-Pearson confidence interval",
    "category": "section",
    "text": "ECPInterval{T<:Real}(\r\n  estimator::Estimator;\r\n  confidence::T = .9,\r\n  k::Integer = 1001\r\n)"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.MMWInterval-Union{Tuple{TE,TR,Function,MathProgBase.SolverInterface.AbstractMathProgSolver}, Tuple{TE}, Tuple{TR}} where TR<:Real where TE<:BinaryTwoStageDesigns.Estimator",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.MMWInterval",
    "category": "method",
    "text": "MMWInterval(\n  estimator::TE,\n  rho0::TR,\n  prior::Function,\n  solver::MathProgBase.AbstractMathProgSolver;\n  confidence::TR = 0.9,\n  ngrid::Integer = 100\n) where {TE<:Estimator,TR<:Real}\n\nExact confidence interval based on ordering induced by estimator minimizing the expected squared width with respect to weight function prior(p::Real).\n\nParameters\n\nParameter Description\nestimator estimator object defining the sample space ordering\nrho0 upper boundary of null hypothesis\nconfidence confidence level of the interval\nsolver MathProgBase solver used for optimization, must support quadratic expressions\nngrid number of equally spaced grid-points on which to check coverage\n\n\n\n"
},

{
    "location": "inference.html#Minimum-mean-width-confidence-interval-1",
    "page": "Inference",
    "title": "Minimum mean width confidence interval",
    "category": "section",
    "text": "MMWInterval(\r\n  estimator::TE,\r\n  rho0::TR,\r\n  prior::Function,\r\n  solver::MathProgBase.AbstractMathProgSolver;\r\n  confidence::TR = 0.9,\r\n  ngrid::Integer = 100\r\n) where {TE<:Estimator,TR<:Real}"
},

]}
