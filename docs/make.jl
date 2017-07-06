using Documenter, BinaryTwoStageDesigns

makedocs(
    format = :html,
    sitename = "BinaryTwoStageDesigns",
    pages = [
        "index.md",
        "Sample Spaces" => "sample_space.md",
        "Design Parameters" => "parameters.md",
        "Two-Stage Designs" => "designs.md",
        "Finding Optimal Designs" => "optimal_designs.md",
        "Inference" => "inference.md"
    ],
    repo = "https://github.com/imbi-heidelberg/BinaryTwoStageDesigns/{commit}"
)
