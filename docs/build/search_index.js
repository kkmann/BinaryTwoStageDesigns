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
    "text": "Depth = 3"
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
    "text": "SampleSpace\n\nAbstract type representing generic sample spaces defining the search space for finding optimal two stage designs.\n\n\n\n"
},

{
    "location": "sample_space.html#BinaryTwoStageDesigns.SimpleSampleSpace-Tuple{Any,T<:Integer}",
    "page": "Sample Spaces",
    "title": "BinaryTwoStageDesigns.SimpleSampleSpace",
    "category": "Method",
    "text": "SimpleSampleSpace{T<:Integer}(n1range, nmax::T; n2min::T = 1, maxnfact::Real = Inf, nmincont::T = 0, maxvariables::T = 500000, GS::Bool = false)\n\nConstructs an object of type SimpleSampleSpace defining the search space for finding optimal two stage designs.\n\nParameters\n\nParameter Description\nn1range possible stage-one sample sizes, can be anything which is convertable to an integer vector\nn2min minimal stage-two sample size\nnmax maximal overall sample size (stage one and two combined)\nmaxnfact nmax must be smaller than maxnfact*n1\nnmincont minimal sample size upon rejection of the null hypothesis\nmaxvariables (approximate) maximal number of binary variables used when finding optimal design; if the sample space is too big this can be used to find approximate solutions to the optimization problem\nGS flag indicating wheather the design should be group-sequential (constant stage two sample size upon continuation)\n\nReturn Value\n\nAn object of type SimpleSampleSpace with the respective parameters.\n\nExamples\n\njulia> SimpleSampleSpace(10:25, 100, n2min = 5)\n\n\n\n"
},

{
    "location": "sample_space.html#Sample-Spaces-1",
    "page": "Sample Spaces",
    "title": "Sample Spaces",
    "category": "section",
    "text": "Sample spaces are used to encode the feasible space of binary two-stage designs for optimization. Currently, only a single option SimpleSampleSpace exists which allows to specify a range of possible stage-one sample sizes, maximal overall sample size and various nicety constraints.SampleSpace\n\nSimpleSampleSpace{T<:Integer}(n1range, nmax::T; n2min::T = 1, maxnfact::Real = Inf, nmincont::T = 0, maxvariables::T = 500000, GS::Bool = false)\n"
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
    "text": "Parameters\n\nAbstract type representing a generic set of paramters for finding optimal two-stage designs.\n\n\n\n"
},

{
    "location": "parameters.html#BinaryTwoStageDesigns.PointAlternative",
    "page": "Design Parameters",
    "title": "BinaryTwoStageDesigns.PointAlternative",
    "category": "Type",
    "text": "PointAlternative <: Parameters\n\nAbstract type representing a generic set of paramters for finding optimal two-stage designs which are characterized by a specific point alternative.\n\n\n\n"
},

{
    "location": "parameters.html#Overview-1",
    "page": "Design Parameters",
    "title": "Overview",
    "category": "section",
    "text": "Objects of type Parameters are used to store parameters required for finding optimal two-stage designs. Every Parameters object contains a sample space which restricts the feasible region of the sample space and implements a score method which is to be optimized.Parameters\n\nPointAlternative"
},

{
    "location": "parameters.html#BinaryTwoStageDesigns.SimpleMinimalExpectedSampleSize",
    "page": "Design Parameters",
    "title": "BinaryTwoStageDesigns.SimpleMinimalExpectedSampleSize",
    "category": "Type",
    "text": "SimpleMinimalExpectedSampleSize{T_samplespace<:SampleSpace} <: PointAlternative\n\nThis type represents a set of parameters for finding optimal two-stage designs minimizing the expected sample size on a point in the parameter space subject to type one and  two error rate constraints.\n\n\n\n"
},

{
    "location": "parameters.html#SimpleMinimalExpectedSampleSize-1",
    "page": "Design Parameters",
    "title": "SimpleMinimalExpectedSampleSize",
    "category": "section",
    "text": "Currently, only a single parameter type (and thus objective function) is implemented. A SimpleMinimalExpectedSampleSize objects holds information about the null hypothesis and a point alternative on which the power is calculated as well as allowing to specify a response rate under which the expected sample size is to be minimized.SimpleMinimalExpectedSampleSize"
},

{
    "location": "designs.html#",
    "page": "Two-Stage Designs",
    "title": "Two-Stage Designs",
    "category": "page",
    "text": ""
},

{
    "location": "designs.html#BinaryTwoStageDesigns.BinaryTwoStageDesign",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.BinaryTwoStageDesign",
    "category": "Type",
    "text": "BinaryTwoStageDesign{T1<:Integer, T2<:Real, PType<:Parameters}\n\nImmutable type representing a binary two-stage design. Two constructors are implemented: BinaryTwoStageDesign{T1<:Integer, T2<:Real}(n::Vector{T1}, c::Vector{T2}) does not require any additional parameters but only vectors n and c of final sample sizes and critical values for x1 = 0:n1, n1 + 1 = length(c) = length(n). The second option BinaryTwoStageDesign{T1<:Integer, T2<:Real, PType<:Parameters}(n::Vector{T1},c::Vector{T2},params::PType) also stores the passed paramter object.\n\nParameters\n\nParameter Description\nn integer vector of length n1 + 1 with sample size n(x1)\nc real vector of length n1 + 1 with critical values c(x1), c(x1) = +/- Inf indicates early stopping and must imply n(x1) = n1\n[params] parameter object\n\nDetails\n\nNote that the test rejects whenever x1 + x2 > c(x1)\n\n\n\n"
},

{
    "location": "designs.html#Base.convert-Tuple{Type{DataFrames.DataFrame},BinaryTwoStageDesigns.BinaryTwoStageDesign}",
    "page": "Two-Stage Designs",
    "title": "Base.convert",
    "category": "Method",
    "text": "convert(::Type{DataFrames.DataFrame}, design::BinaryTwoStageDesign)\n\nConvert a design to a DataFrame with columns x1, n, c.\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.interimsamplesize-Tuple{BinaryTwoStageDesigns.BinaryTwoStageDesign}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.interimsamplesize",
    "category": "Method",
    "text": "interimsamplesize(design::BinaryTwoStageDesign)\n\nReturn inter sample size n1 of a design.\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.parameters-Tuple{BinaryTwoStageDesigns.BinaryTwoStageDesign}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.parameters",
    "category": "Method",
    "text": "parameters(design::BinaryTwoStageDesign)\n\nReturn stored parameter object of a design.\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.samplesize-Tuple{BinaryTwoStageDesigns.BinaryTwoStageDesign}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.samplesize",
    "category": "Method",
    "text": "samplesize(design::BinaryTwoStageDesign)\n\nReturn vector of final sample sizes for x1 = 0:n1, or via samplesize{T<:Integer}(design::BinaryTwoStageDesign, x1::T) only for specific x1.\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.criticalvalue-Tuple{BinaryTwoStageDesigns.BinaryTwoStageDesign}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.criticalvalue",
    "category": "Method",
    "text": "criticalvalue(design::BinaryTwoStageDesign)\n\nReturn vector of critical values for x1 = 0:n1, or via criticalvalue(design::BinaryTwoStageDesign, x1) only for specific x1.\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.power-Tuple{BinaryTwoStageDesigns.BinaryTwoStageDesign,T<:Real}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.power",
    "category": "Method",
    "text": "power{T<:Real}(design::BinaryTwoStageDesign, p::T)\n\nCompute the power of a given design.\n\nParameters\n\nParameter Description\ndesign a BinaryTwoStageDesign\np response probability\n\nDetails\n\nA conditional power method is also implemented as power{T1<:Integer, T2<:Real}(design::BinaryTwoStageDesign, x1::T1, p::T2) where x1 is the stage-one number of responses.\n\nReturn Value\n\nConditional power of design at p given x1 responses in stage one.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> params = SimpleMinimalExpectedSampleSize(ss, .2, .4, .05, .2, .4)\njulia> design, res = getoptimaldesign(params, solver = Gurobi.GurobiSolver())\njulia> power(design, .4)\njulia> power(design, 0, .4)\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.stoppingforfutility-Tuple{BinaryTwoStageDesigns.BinaryTwoStageDesign,T<:Real}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.stoppingforfutility",
    "category": "Method",
    "text": "stoppingforfutility{T<:Real}(design::BinaryTwoStageDesign, p::T)\n\nCompute probability of stopping early for futility of a given design.\n\nParameters\n\nParameter Description\ndesign a BinaryTwoStageDesign\np response probability\n\nReturn Value\n\nProbablity of stopping-for-futility\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> params = SimpleMinimalExpectedSampleSize(ss, .2, .4, .05, .2, .4)\njulia> design, res = getoptimaldesign(params, solver = Gurobi.GurobiSolver())\njulia> stoppingforfutility(design, .2)\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.test-Tuple{BinaryTwoStageDesigns.BinaryTwoStageDesign,T<:Integer,T<:Integer}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.test",
    "category": "Method",
    "text": "test{T<:Integer}(design::BinaryTwoStageDesign, x1::T, x2::T)::Bool\n\nBinary test decision of design when observing x1 responses in stage one and x2 responses in stage two.\n\nParameters\n\nParameter Description\ndesign a BinaryTwoStageDesign\nx1 number of stage-one responses\nx2 number of stage-two responses\n\nReturn Value\n\nBoolean, true if x1 + x2 > c(x1), i.e., if the null hypothesis can be rejected.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> params = SimpleMinimalExpectedSampleSize(ss, .2, .4, .05, .2, .4)\njulia> design, res = getoptimaldesign(params, solver = Gurobi.GurobiSolver())\njulia> test(design, 0, 0)\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.simulate-Tuple{BinaryTwoStageDesigns.BinaryTwoStageDesign,T2<:Real,T1<:Integer}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.simulate",
    "category": "Method",
    "text": "simulate{T1<:Integer, T2<:Real}(design::BinaryTwoStageDesign, p::T2, nsim::T1)\n\nSimulate nsim trials of design design with true response probability p.\n\nParameters\n\nParameter Description\ndesign a BinaryTwoStageDesign\np true response probability\nnsim number of simulated trials\n\nReturn Value\n\nA DataFrame with columns x1 (stage one responses) n (overall sample size) c (critical value), x2 (stage-two responses), and rejectedH0 (test decision).\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> params = SimpleMinimalExpectedSampleSize(ss, .2, .4, .05, .2, .4)\njulia> design, res = getoptimaldesign(params, solver = Gurobi.GurobiSolver())\njulia> df = simulate(design, .3, 1000)\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.score-Tuple{BinaryTwoStageDesigns.BinaryTwoStageDesign,BinaryTwoStageDesigns.Parameters}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.score",
    "category": "Method",
    "text": "score(design::BinaryTwoStageDesign, params::Parameters)::Real\n\nEvaluates the score of a design for a given set of parameters.\n\nParameters\n\nParameter Description\ndesign a BinaryTwoStageDesign\nparams a parameters object holding the required information about the score\n\nDetails\n\nA method score(design::BinaryTwoStageDesign)::Real is also implemented which uses the parameter object stored within the design object after optimization. Note that this is only possible if the design was created as result of calling getoptimaldesign().\n\nReturn Value\n\nA single real value, i.e., the score.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> params = SimpleMinimalExpectedSampleSize(ss, .2, .4, .05, .2, .4)\njulia> design, res = getoptimaldesign(params, solver = Gurobi.GurobiSolver())\njulia> score(design, params)\njulia> score(design)\n\n\n\n"
},

{
    "location": "designs.html#BinaryTwoStageDesigns.jeffreysprior-Tuple{BinaryTwoStageDesigns.BinaryTwoStageDesign}",
    "page": "Two-Stage Designs",
    "title": "BinaryTwoStageDesigns.jeffreysprior",
    "category": "Method",
    "text": "jeffreysprior(design::BinaryTwoStageDesign)\n\nComputes the Jeffreys prior of any given design.\n\nParameters\n\nParameter Description\ndesign a BinaryTwoStageDesign\n\nReturn Value\n\nA function with signature prior{T<:Real}(p::T)::Real where p is the response probability and prior(p) the PDF of the Jeffres prior at p.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> params = SimpleMinimalExpectedSampleSize(ss, .2, .4, .05, .2, .4)\njulia> design, res = getoptimaldesign(params, solver = Gurobi.GurobiSolver())\njulia> f = jeffreysprior(design)\njulia> f(.5)\n\n\n\n"
},

{
    "location": "designs.html#Distributions.pdf-Tuple{BinaryTwoStageDesigns.BinaryTwoStageDesign,T1<:Integer,T1<:Integer,T2<:Real}",
    "page": "Two-Stage Designs",
    "title": "Distributions.pdf",
    "category": "Method",
    "text": "pdf{T1<:Integer, T2<:Real}(design::BinaryTwoStageDesign, x1::T1, x2::T1, p::T2)\n\nProbability density function of (x1, x2) responses in stage one and two respectively under design given response rate p.\n\nParameters\n\nParameter Description\ndesign a BinaryTwoStageDesign\nx1 number of stage-one responses\nx2 number of stage-two responses\np response probability\n\nReturn Value\n\nPDF of (x1, x2) given design, p\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> params = SimpleMinimalExpectedSampleSize(ss, .2, .4, .05, .2, .4)\njulia> design, res = getoptimaldesign(params, solver = Gurobi.GurobiSolver())\njulia> pdf(design, 0, 0, .4)\n\n\n\n"
},

{
    "location": "designs.html#Binary-Two-Stage-Designs-1",
    "page": "Two-Stage Designs",
    "title": "Binary Two-Stage Designs",
    "category": "section",
    "text": "BinaryTwoStageDesign{T1<:Integer, T2<:Real, PType<:Parameters}\n\nconvert(::Type{DataFrames.DataFrame}, design::BinaryTwoStageDesign)\n\ninterimsamplesize(design::BinaryTwoStageDesign)\n\nparameters(design::BinaryTwoStageDesign)\n\nsamplesize(design::BinaryTwoStageDesign)\n\ncriticalvalue(design::BinaryTwoStageDesign)\n\npower{T<:Real}(design::BinaryTwoStageDesign, p::T)\n\nstoppingforfutility{T<:Real}(design::BinaryTwoStageDesign, p::T)\n\ntest{T<:Integer}(design::BinaryTwoStageDesign, x1::T, x2::T)\n\nsimulate{T1<:Integer, T2<:Real}(design::BinaryTwoStageDesign, p::T2, nsim::T1)\n\nscore(design::BinaryTwoStageDesign, params::Parameters)\n\njeffreysprior(design::BinaryTwoStageDesign)\n\npdf{T1<:Integer, T2<:Real}(design::BinaryTwoStageDesign, x1::T1, x2::T1, p::T2)"
},

{
    "location": "optimal_designs.html#",
    "page": "Finding Optimal Designs",
    "title": "Finding Optimal Designs",
    "category": "page",
    "text": ""
},

{
    "location": "optimal_designs.html#BinaryTwoStageDesigns.getoptimaldesign-Tuple{T<:Integer,BinaryTwoStageDesigns.Parameters,TS<:MathProgBase.SolverInterface.AbstractMathProgSolver}",
    "page": "Finding Optimal Designs",
    "title": "BinaryTwoStageDesigns.getoptimaldesign",
    "category": "Method",
    "text": "getoptimaldesign{T<:Integer, TS<:MathProgBase.AbstractMathProgSolver}(\n    n1::T,\n    parameters::Parameters,\n    solver::TS;\n    VERBOSE::Integer = 0\n)\n\nFind the optimal two-stage design for given parameters and fixed stage-one sample size.\n\nParameters\n\nParameter Description\nn1 stage-one sample size\nnparameters paramters object defining the optimality criterion\nsolver MathProgBase solver used for optimization\nVERBOSE control verbosity during optimization\n\nReturn Value\n\nBinaryTwoStageDesign object optimized for given parameter set.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())\n\n\n\n"
},

{
    "location": "optimal_designs.html#BinaryTwoStageDesigns.getoptimaldesign-Tuple{BinaryTwoStageDesigns.Parameters,TS<:MathProgBase.SolverInterface.AbstractMathProgSolver}",
    "page": "Finding Optimal Designs",
    "title": "BinaryTwoStageDesigns.getoptimaldesign",
    "category": "Method",
    "text": "getoptimaldesign{TS<:MathProgBase.AbstractMathProgSolver}(\n    parameters::Parameters,\n    solver::TS;\n    VERBOSE::Integer = 1\n)\n\nFind the optimal two-stage design for given parameters (optimizes over n1 as well).\n\nParameters\n\nParameter Description\nnparameters paramters object defining the optimality criterion\nsolver MathProgBase solver used for optimization\nVERBOSE control verbosity during optimization\n\nReturn Value\n\nTuple, first element is the BinaryTwoStageDesign object optimized for given parameter set and the second element is a dictionary with n1, scores and respective optimal designs.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design, res = getoptimaldesign(params, solver = Gurobi.GurobiSolver())\n\n\n\n"
},

{
    "location": "optimal_designs.html#Optimal-Two-Stage-Designs-1",
    "page": "Finding Optimal Designs",
    "title": "Optimal Two-Stage Designs",
    "category": "section",
    "text": "The principle technical background in available at arxiv.org.getoptimaldesign{T<:Integer, TS<:MathProgBase.AbstractMathProgSolver}(\n    n1::T,\n    parameters::Parameters,\n    solver::TS;\n    VERBOSE::Integer = 0\n)\n\ngetoptimaldesign{T<:Integer, TS<:MathProgBase.AbstractMathProgSolver}(\n    parameters::Parameters,\n    solver::TS;\n    VERBOSE::Integer = 1\n)"
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
    "text": "Pages = [\"inference.md\"]\nDepth = 3"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.BinaryTwoStageDesignEstimator",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.BinaryTwoStageDesignEstimator",
    "category": "Type",
    "text": "BinaryTwoStageDesignEstimator\n\nAbstract base type for all estimators.\n\n\n\n"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.estimate-Tuple{BinaryTwoStageDesigns.BinaryTwoStageDesignEstimator,T<:Integer,T<:Integer}",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.estimate",
    "category": "Method",
    "text": "estimate{T<:Integer}(estimator::BinaryTwoStageDesignEstimator, x1::T, x2::T)\n\nEstimate the response rate from observed x1 and x2.\n\nParameters\n\nParameter Description\nestimator any BinaryTwoStageDeisgnEstimator object\nx1 stage-one responses\nx2 stage-two responses\n\nReturn Value\n\nReal, estimated response rate.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())\njulia> est = MLE(design)\njulia> estimate(est, 0, 0)\n\n\n\n"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.bias-Tuple{BinaryTwoStageDesigns.BinaryTwoStageDesignEstimator,T<:Real}",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.bias",
    "category": "Method",
    "text": "bias{T<:Real}(estimator::BinaryTwoStageDesignEstimator, p::T)\n\nBias of estimator given response rate p.\n\nParameters\n\nParameter Description\nestimator any BinaryTwoStageDeisgnEstimator object\np0 upper boundary of null hypothesis\n\nReturn Value\n\nReal, bias given p.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())\njulia> est = MLE(design)\njulia> bias(est, .3)\n\n\n\n"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.rmse-Tuple{BinaryTwoStageDesigns.BinaryTwoStageDesignEstimator,T<:Real}",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.rmse",
    "category": "Method",
    "text": "rmse{T<:Real}(estimator::BinaryTwoStageDesignEstimator, p::T)\n\nRoot mean squared error of estimator given response rate p.\n\nParameters\n\nParameter Description\nestimator any BinaryTwoStageDeisgnEstimator object\np0 upper boundary of null hypothesis\n\nReturn Value\n\nReal, RMSE given p.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())\njulia> est = MLE(design)\njulia> rmse(est, .3)\n\n\n\n"
},

{
    "location": "inference.html#Point-estimation-1",
    "page": "Inference",
    "title": "Point estimation",
    "category": "section",
    "text": "BinaryTwoStageDesignEstimator\n\nestimate{T<:Integer}(estimator::BinaryTwoStageDesignEstimator, x1::T, x2::T)\n\nbias{T<:Real}(estimator::BinaryTwoStageDesignEstimator, p::T)\n\nrmse{T<:Real}(estimator::BinaryTwoStageDesignEstimator, p::T)"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.MaximumLikelihoodEstimator",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.MaximumLikelihoodEstimator",
    "category": "Type",
    "text": "MaximumLikelihoodEstimator\n\nMaximumLikelihoodEstimator(design::BinaryTwoStageDesign)\n\nSimple maximum likelihood estimator for response rate p.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())\njulia> est = MaximumLikelihoodEstimator(design)\n\n\n\n"
},

{
    "location": "inference.html#Maximum-likelihood-estimator-1",
    "page": "Inference",
    "title": "Maximum likelihood estimator",
    "category": "section",
    "text": "MaximumLikelihoodEstimator"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.RaoBlackwellizedEstimator",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.RaoBlackwellizedEstimator",
    "category": "Type",
    "text": "RaoBlackwellizedEstimator <: BinaryTwoStageDesignEstimator\n\nRaoBlackwellizedEstimator(design::BinaryTwoStageDesign)\n\nUnbiased estimator for response rate p see also:\n\nKunzmann K, Kieser M. Point estimation and p‐values in phase II adaptive two‐stage designs with a binary endpoint. Statistics in medicine. 2017 Mar 15;36(6):971-84.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())\njulia> est = RaoBlackwellizedEstimator(design)\n\n\n\n"
},

{
    "location": "inference.html#Unbiased-estimator-1",
    "page": "Inference",
    "title": "Unbiased estimator",
    "category": "section",
    "text": "RaoBlackwellizedEstimator"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.CompatibleEstimator",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.CompatibleEstimator",
    "category": "Type",
    "text": "CompatibleEstimator <: BinaryTwoStageDesignEstimator\n\nCompatibleEstimator{TS<:MathProgBase.AbstractMathProgSolver}(\n    design::BinaryTwoStageDesign,\n    solver::TS;\n    prior::Function = jeffreysprior(design),\n    k = 100\n)\n\nCompatible estimator minimizing expected MSE for response rate p see also:\n\nKunzmann K, Kieser M. Point estimation and p‐values in phase II adaptive two‐stage designs with a binary endpoint. Statistics in medicine. 2017 Mar 15;36(6):971-84.\n\nParameters\n\nParameter Description\ndesign BinaryTwoStageDesign\nsolver MathProgBase solver, must support quadratic expressions\nprior weight function for MSE values at different p, must be of form f(p::Real)::Real\nk number of equally spaced grid-points for evaluation of MSE and prior\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, Gurobi.GurobiSolver())\njulia> est = CompatibleEstimator(design, Gurobi.GurobiSolver())\n\n\n\n"
},

{
    "location": "inference.html#Optimal-compatible-estimator-1",
    "page": "Inference",
    "title": "Optimal compatible estimator",
    "category": "section",
    "text": "CompatibleEstimator"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.p-Tuple{BinaryTwoStageDesigns.BinaryTwoStageDesignEstimator,T1<:Integer,T1<:Integer,T2<:Real}",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.p",
    "category": "Method",
    "text": "p{T1<:Integer, T2<:Real}(estimator::BinaryTwoStageDesignEstimator, x1::T1, x2::T1, p0::T2)\n\nCompute the p value after observing (x1, x2) for null hypothesis H0: p <= p0 with respect to ordering induced by estimator.\n\nParameters\n\nParameter Description\nestimator any BinaryTwoStageDeisgnEstimator object\nx1 stage-one responses\nx2 stage-two responses\np0 upper boundary of null hypothesis\n\nReturn Value\n\nReal, p value.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())\njulia> est = MLE(design)\njulia> p(est, 0, 0, .2)\n\n\n\n"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.incompatibleoutcomes-Tuple{BinaryTwoStageDesigns.BinaryTwoStageDesignEstimator}",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.incompatibleoutcomes",
    "category": "Method",
    "text": "incompatibleoutcomes(estimator::BinaryTwoStageDesignEstimator)\n\nFind outcomes where the induced p value implies different decision than the underlying design\n\nParameters\n\nParameter Description\nestimator any BinaryTwoStageDeisgnEstimator object\n\nReturn Value\n\nArray with respective outcomes\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())\njulia> est = MLE(design)\njulia> incompatibleoutcomes(est)\n\n\n\n"
},

{
    "location": "inference.html#P-values-1",
    "page": "Inference",
    "title": "P values",
    "category": "section",
    "text": "p{T1<:Integer, T2<:Real}(estimator::BinaryTwoStageDesignEstimator, x1::T1, x2::T1, p0::T2)\n\nincompatibleoutcomes(estimator::BinaryTwoStageDesignEstimator)"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.ConfidenceInterval",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.ConfidenceInterval",
    "category": "Type",
    "text": "ConfidenceInterval\n\nAbstract base type for all confidence interval types.\n\n\n\n"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.limits",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.limits",
    "category": "Function",
    "text": "limits{T<:Integer}(ci::ConfidenceInterval, x1::T, x2::T)\n\nReturn the confidence interval limits for observed x1 and x2.\n\nParameters\n\nParameter Description\nci any ConfidenceInterval object\nx1 stage-one responses\nx2 stage-two responses\n\nReturn Value\n\nTwo element Real vector of limits [lower, upper].\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())\njulia> est = MLE(design)\njulia> ci = ClopperPearsonConfidenceInterval(est, confidence = .9)\njulia> limits(ci, 0, 0)\n\n\n\n"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.coverage",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.coverage",
    "category": "Function",
    "text": "coverage{T<:Real}(ci::ConfidenceInterval, p::T; orientation = \"overall\")\n\nReturn coverage of given confidence interval and response rate p.\n\nParameters\n\nParameter Description\nci any ConfidenceInterval object\np response rate\norientation string indicating the coverage type - \"overall\", \"lower\", or \"upper\"\n\nReturn Value\n\nCoverage probability given p.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())\njulia> est = MLE(design)\njulia> ci = ClopperPearsonConfidenceInterval(est, confidence = .9)\njulia> coverage(ci, .5)\n\n\n\n"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.meanwidth",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.meanwidth",
    "category": "Function",
    "text": "meanwidth{T<:Real}(ci::ConfidenceInterval, p::T)\n\nMean width of given confidence interval and response rate p.\n\nParameters\n\nParameter Description\nci any ConfidenceInterval object\np response rate\n\nReturn Value\n\nMean width given p.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())\njulia> est = MLE(design)\njulia> ci = ClopperPearsonConfidenceInterval(est, confidence = .9)\njulia> meanwidth(ci, .5)\n\n\n\n"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.meaninterval",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.meaninterval",
    "category": "Function",
    "text": "meaninterval{T<:Real}(ci::ConfidenceInterval, p::T)\n\nMean interval (average limits) of given confidence interval and response rate p.\n\nParameters\n\nParameter Description\nci any ConfidenceInterval object\np response rate\n\nReturn Value\n\nMean interval given p.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())\njulia> est = MLE(design)\njulia> ci = ClopperPearsonConfidenceInterval(est, confidence = .9)\njulia> meaninterval(ci, .5)\n\n\n\n"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.findinconsistencies",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.findinconsistencies",
    "category": "Function",
    "text": "findinconsistencies{T<:Real}(ci::ConfidenceInterval)\n\nReturn outcomes where the confidence interval contradicts the designs test decision.\n\nParameters\n\nParameter Description\nci any ConfidenceInterval object\np response rate\n\nReturn Value\n\nMean interval given p.\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, solver = Gurobi.GurobiSolver())\njulia> est = MLE(design)\njulia> ci = ClopperPearsonConfidenceInterval(est, confidence = .9)\njulia> findinconsistencies(ci)\n\n\n\n"
},

{
    "location": "inference.html#Confidence-intervals-1",
    "page": "Inference",
    "title": "Confidence intervals",
    "category": "section",
    "text": "ConfidenceInterval\n\nlimits\n\ncoverage\n\nmeanwidth\n\nmeaninterval\n\nfindinconsistencies"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.NaiveClopperPearsonConfidenceInterval",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.NaiveClopperPearsonConfidenceInterval",
    "category": "Type",
    "text": "NaiveClopperPearsonConfidenceInterval <: ConfidenceInterval\n\nNaiveClopperPearsonConfidenceInterval{T<:Real}(\n    design::BinaryTwoStageDesign;\n    confidence::T = .9\n)\n\nNaive Clopper-Pearson confidence interval using default ordering.\n\nParameters\n\nParameter Description\nestimator estimator object defining the sample space ordering\nconfidence confidence level of the interval\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, Gurobi.GurobiSolver())\njulia> est = MaximumLikelihoodEstimator(design, Gurobi.GurobiSolver())\njulia> ci = NaiveClopperPearsonConfidenceInterval(est, confidence = .9)\n\n\n\n"
},

{
    "location": "inference.html#Naive-Clopper-Pearson-confidence-interval-1",
    "page": "Inference",
    "title": "Naive Clopper-Pearson confidence interval",
    "category": "section",
    "text": "NaiveClopperPearsonConfidenceInterval"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.ClopperPearsonConfidenceInterval",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.ClopperPearsonConfidenceInterval",
    "category": "Type",
    "text": "ClopperPearsonConfidenceInterval <: ConfidenceInterval\n\nClopperPearsonConfidenceInterval{T<:Real}(\n    estimator::BinaryTwoStageDesignEstimator;\n    confidence::T = .9,\n    k::Integer = 1001\n)\n\nExact Clopper-Pearson type confidence interval based on ordering induced by estimator.\n\nParameters\n\nParameter Description\nestimator estimator object defining the sample space ordering\nconfidence confidence level of the interval\nk number of equally spaced grid-points for invertign the test\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, Gurobi.GurobiSolver())\njulia> est = MaximumLikelihoodEstimator(design, Gurobi.GurobiSolver())\njulia> ci = ClopperPearsonConfidenceInterval(est, confidence = .9)\n\n\n\n"
},

{
    "location": "inference.html#Clopper-Pearson-confidence-interval-1",
    "page": "Inference",
    "title": "Clopper-Pearson confidence interval",
    "category": "section",
    "text": "ClopperPearsonConfidenceInterval"
},

{
    "location": "inference.html#BinaryTwoStageDesigns.MinimumMeanWidthConfidenceInterval",
    "page": "Inference",
    "title": "BinaryTwoStageDesigns.MinimumMeanWidthConfidenceInterval",
    "category": "Type",
    "text": "MinimumMeanWidthConfidenceInterval <: ConfidenceInterval\n\nMinimumMeanWidthConfidenceInterval{TS<:MathProgBase.AbstractMathProgSolver}(\n    estimator::BinaryTwoStageDesignEstimator,\n    rho0::Float64,\n    prior::Function,\n    solver::TS;\n    confidence::Float64 = .9,\n    ngrid::Int64 = 100\n)\n\nExact confidence interval based on ordering induced by estimator minimizing the expected squared width with respect to weight function prior(p::Real).\n\nParameters\n\nParameter Description\nestimator estimator object defining the sample space ordering\nrho0 upper boundary of null hypothesis\nconfidence confidence level of the interval\nsolver MathProgBase solver used for optimization, must support quadratic expressions\nngrid number of equally spaced grid-points on which to check coverage\n\nExamples\n\njulia> ss = SimpleSampleSpace(10:25, 100, n2min = 5)\njulia> interimsamplesize(ss)\njulia> design = getoptimaldesign(15, params, Gurobi.GurobiSolver())\njulia> est = MaximumLikelihoodEstimator(design, Gurobi.GurobiSolver())\njulia> ci = ClopperPearsonConfidenceInterval(est, confidence = .9)\n\n\n\n"
},

{
    "location": "inference.html#Minimum-mean-width-confidence-interval-1",
    "page": "Inference",
    "title": "Minimum mean width confidence interval",
    "category": "section",
    "text": "MinimumMeanWidthConfidenceInterval"
},

]}
