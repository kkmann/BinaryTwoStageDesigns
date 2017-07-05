using Documenter, BinaryTwoStageDesigns

makedocs(
    format = :html,
    sitename = "BinaryTwoStageDesigns",
    pages = [
        "index.md",
        "Sample Spaces" => "sample_space.md",
        "Design Parameters" => "parameters.md",
        "BinaryTwoStageDesign" => "designs.md"
    ],
    repo = "https://github.com/imbi-heidelberg/BinaryTwoStageDesigns/{commit}"
    )
