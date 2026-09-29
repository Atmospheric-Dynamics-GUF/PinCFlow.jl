# Get the version number of the old release.
old_version = try
    VersionNumber(read(`git describe --tags --abbrev=0`, String))
catch
    error("No previous release was found!")
end
(; major, minor, patch) = old_version

# Get the changes since the old release.
run(`git fetch origin main`)
changes = split(
    read(`git log v$old_version..origin/main --pretty=format:"  - %s"`, String),
    "\n",
)
commit_count = length(changes)
commit_count == 0 &&
    error("No commits have been added since the previous release!")

# Extract the PR numbers.
pr_numbers = Set{Int64}()
for change in changes
    numbers = Tuple(
        parse(Int64, m.captures[1]) for m in eachmatch(r"\(#(\d+)\)", change)
    )
    length(numbers) == 0 &&
        println("WARNING: Commit message without PR reference detected!")
    push!(pr_numbers, numbers...)
end

# Extract the PR labels.
pr_labels = Set{String}()
for pr_number in pr_numbers
    labels = read(
        `gh pr view $pr_number --json labels --jq ".labels[].name"`,
        String,
    )
    labels == "" && error("PR #$pr_number does not have a label!")
    push!(pr_labels, split(labels, "\n")...)
end

# Determine the version number of the new release.
new_version = if "breaking change" in pr_labels
    VersionNumber(major + 1, 0, 0)
elseif "feature change" in pr_labels
    VersionNumber(major, minor + 1, 0)
elseif "bug fix" in pr_labels
    VersionNumber(major, minor, patch + 1)
else
    error("No new release required!")
end

# Format the changes.
changes = replace(
    join(changes, "\n\n"),
    r"\(#(\d+)\)" =>
        s"([#\1](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/\1))",
)


# Update the Project.toml files.
for (pattern, file) in zip(
    ("PinCFlow", "PinCFlow", "version"),
    ("docs/Project.toml", "examples/Project.toml", "Project.toml"),
)
    write(
        file,
        replace(
            read(file, String),
            Regex("$pattern = ") *
            r"\"\d+\.\d+\.\d+\"" => "$pattern = \"$new_version\"",
        ),
    )
end

# Update the changelog.
old_header = "## Release $old_version"
new_header = "## Release $new_version"
changelog = read("NEWS.md", String)
if !occursin(new_header, changelog)
    write(
        "NEWS.md",
        replace(
            changelog,
            Regex(old_header) => "$new_header\n\n$changes\n\n$old_header",
        ),
    )
end
