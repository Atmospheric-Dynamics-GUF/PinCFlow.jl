"""This Fortran code formatter unifies spacing, subunit indents and line \
    breaks in Fortran files. The spacing in comments and strings is not \
    changed. Comments between continuation lines are moved below the last of \
    these lines. Comments between quotes (i.e. in Fortran strings that contain \
    line breaks) are removed (a warning is given in such cases). The output \
    does not contain any extraneous whitespaces or consecutive blank lines \
    (outside of quotes and comments). The user may specify input and output \
    directories, as well as the lengths of the spacing and indent strings."""

import os
import re
import sys

# Set directories.
if len(sys.argv) >= 2:
  input_directory = sys.argv[1]
else:
  input_directory = "."
if len(sys.argv) >= 3:
  output_directory = sys.argv[2]
else:
  output_directory = "."

# Create output directory if not existing.
if not os.path.exists(output_directory):
  os.makedirs(output_directory)

# Set parameters.
if len(sys.argv) >= 4:
  subunit_indent = int(sys.argv[3])
else:
  subunit_indent = 2
if len(sys.argv) == 5:
  linebreak_indent = int(sys.argv[4])
else:
  linebreak_indent = 4

# More than two parameters are not allowed.
if len(sys.argv) > 5:
  sys.exit("Too many input parameters!")

# Define indent strings.
subunit_indent_string = subunit_indent * " "
linebreak_indent_string = linebreak_indent * " "

# Collect operators.
operators = re.compile(r"(\+|\-|\*|\/|\=|\>|\<|\.eq\.|\.ne\.|\.gt\.|\.lt\." \
    r"|\.ge\.|\.le\.|\.and\.|\.or\.|\.not\.|\.eqv\.|\.neqv\.|\.true\." \
    r"|\.false\.)", flags = re.IGNORECASE)

# Collect indent triggers.
positive_indent_triggers = re.compile("".join((r"^(\d+\s*)?(\w+\s*:\s*)?(", \
    r"|".join((r"if\s*\(.*\)\s*then", \
    r"else(\s*if\s*\(.*\)\s*then)?", \
    r"do(\s+\w+.*)?", \
    r"select\s*(case|rank|type)\s*\(.*\)", \
    r"((case|rank|type\s+is|class\s+is)\s*(\(.*\)|default)|class\s+default)", \
    r"associate\s*\(.*\)", \
    r"block", \
    r"(?!end)(\w+\s+)*\s*subroutine\s+\w+\s*(\(.*\))?", \
    r"(?!end)(\w+\s+)*\s*function\s+\w+\s*(\(.*\))?(\s*result\s*\(\w+\))?", \
    r"module\s+\w+", \
    r"submodule\s*\(\w+\)\s*\w+", \
    r"type(\s*,\s*(bind\s*\(.*\)|extends\s*\(.*\)|abstract|public|private))*" \
    r"(\s*::\s*|\s+)\w+", \
    r"program\s+\w+", \
    r"(?!end)(\w+\s+)*\s*interface(\s+\w+|\s+(operator|assignment)" \
    r"\s*\(.*\))?", \
    r"enum(\s*,\s*(bind\s*\(.*\)))?((\s*::\s*|\s+)\w+)?")), r")$")), \
    flags = re.IGNORECASE)
negative_indent_triggers = re.compile("".join((r"^(\d+\s*)?((", \
    r"|".join((r"else(\s*if\s*\(.*\)\s*then)?", \
    r"end\s*if", \
    r"end\s*do", \
    r"((case|rank|type\s+is|class\s+is)\s*(\(.*\)|default)|" \
    r"class\s+default)", \
    r"end\s*select", \
    r"end\s*associate", \
    r"end\s*block", \
    r"end\s*subroutine", \
    r"end\s*function", \
    r"end\s*module", \
    r"end\s*submodule", \
    r"end\s*type", \
    r"end\s*program", \
    r"end\s*interface(\s+\w+|\s+(operator|assignment)\s*\(.*\))?", \
    r"end\s*enum")), r")(\s+\w+)?|end)$")), flags = re.IGNORECASE)

# Define expressions for line break detection.
linebreak_comments = re.compile(r"(\n[^!\n]+& *(!.*)?\n)" \
    r"\n*((!.*\n)+)\n*(&?[^!\n]+&? *(!.*)?\n)")
linebreak_code = re.compile(r"(\n[^!\n]+)& *((!.*)?)" \
    r"\n&?([^!\n]+&?) *((!.*)?\n)")

# Define expressions for line break points.
operator_linebreak_point = re.compile(r"^(\s*)(\S.*\s+)((\+|\-|\*|\/" \
    r"|\=|\<|\>)+\s+[^\s\+\-\*\/\=\<\>]+\s*)$")
whitespace_linebreak_point = re.compile(r"^(\s*)(\S.*\s+)(\S+\s*)$")

# Iterate over files in the input directory.
for file in os.listdir(input_directory):
  if re.search(r"\.(f|for|ftn|f90|f95|f03|fpp)$", file, flags = re.IGNORECASE) \
      is None:
    continue

  # Read the input file.
  with open("/".join((input_directory, file)), "r") as input:
    code = input.read()

  # Replace tabs with whitespaces.
  code = code.expandtabs()

  # Remove leading and trailing whitespaces.
  code = "\n".join((entry.strip(" ") for entry in code.split("\n")))

  # Remove consecutive blank lines.
  code = re.sub(r"\n\n+", r"\n\n", code)

  # Extract everything between quotes.
  comment = False
  single_quotes = False
  double_quotes = False
  quotes = []
  modified_code = ""
  for character in code:
    modified_code = "".join((modified_code, character))
    if character == "\n":
      comment = False
    elif not (single_quotes or double_quotes) and character == "!":
      comment = True
    if not (comment or double_quotes) and character == "'":
      single_quotes = not single_quotes
      if single_quotes:
        quotes.append("")
    elif not (comment or single_quotes) and character == '"':
      double_quotes = not double_quotes
      if double_quotes:
        quotes.append("")
    elif single_quotes or double_quotes:
      quotes[- 1] = "".join((quotes[- 1], character))
      modified_code = modified_code[:(- 1)]
  code = modified_code

  # Join continuation lines between quotes. Comments between continuation
  # lines are removed.
  for (index, entry) in enumerate(quotes):
    if re.search(r"\n( *!.*\n)+", entry):
      entry = re.sub(r"\n( *!.*\n)+", r"\n", entry)
      print("WARNING: Comments between quotes have been removed!")
    entry = re.sub(r"( *)& *\n *&( *)", r"\1\2", entry)
    entry = re.sub(r"( *)& *\n *", r"\1", entry)
    quotes[index] = entry

  # Join continuation lines. Comments in continuation lines are moved into the
  # first line. Comments between continuation lines are moved below the last
  # line.
  while linebreak_comments.search(code):
    code = linebreak_comments.sub(r"\1\5\3", code)
  while linebreak_code.search(code):
    code = linebreak_code.sub(r"\1\4 \2 \5", code)

  # Separate code and comments.
  comments = ["".join(("!", line.split("!", 1)[1])) if len(line.split("!", 1)) \
      == 2 else "" for line in code.split("\n")]
  code = "\n".join([line.split("!", 1)[0] for line in code.split("\n")])

  # Adjust whitespaces directly outside of parantheses.
  code = re.sub(r" *\(", "(", code)
  code = re.sub(r"\) *", ") ", code)
  code = re.sub(r" *\[", "[", code)
  code = re.sub(r"\] *", "] ", code)

  # Remove whitespaces around percentage sign.
  code = re.sub(r" *\% *", "%", code)

  # Add whitespaces around single operators.
  code = operators.sub(r" \1 ", code)

  # Adjust whitespaces in power notation.
  code = re.sub(r"\b([\d.]+) *(e|d) *((\+|\-)?) *([\d.]+)\b", r"\1\2\3\5", \
      code, flags = re.IGNORECASE)

  # Adjust whitespaces at combined operators.
  code = re.sub(r"\* +\*", "**", code)
  code = re.sub(r"\= +\=", "==", code)
  code = re.sub(r"\/ +\=", "/=", code)
  code = re.sub(r"\< +\=", "<=", code)
  code = re.sub(r"\> +\=", ">=", code)
  code = re.sub(r"\= +\>", "=>", code)
  code = re.sub(r"\/ +\/", "//", code)

  # Adjust whitespaces around colons.
  code = re.sub(r" *: *", ":", code)
  code = re.sub(r"(\n(\d+ *)?\w+:)(\w+)", r"\1 \3", code)
  code = "::".join((" ".join(("", entry, "")) if entry.count("(") \
      == entry.count(")") else entry for entry in code.split("::")))

  # Adjust whitespaces around commas and semicolons.
  code = re.sub(r" *, *", ", ", code)
  code = re.sub(r" *; *", "; ", code)

  # Remove whitespaces directly inside parantheses.
  code = re.sub(r"\( +", "(", code)
  code = re.sub(r" +\)", ")", code)
  code = re.sub(r"\(\/ +", "(/", code)
  code = re.sub(r" +\/\)", "/)", code)
  code = re.sub(r"\[ +", "[", code)
  code = re.sub(r" +\]", "]", code)

  # Remove extraneous whitespaces.
  code = re.sub(r" +", " ", code)

  # Put code and comments back together.
  code = code.split("\n")
  code = "\n".join((" ".join((code[index].rstrip(" "), comments[index])) if \
      comments[index] and code[index].rstrip(" ") else "".join((code[index], \
      comments[index])) for index in range(len(code))))

  # Remove leading and trailing whitespaces.
  code = "\n".join((entry.strip(" ") for entry in code.split("\n")))

  # Add indents.
  code = code.split("\n")
  indent = ""
  for (index, line) in enumerate(code):
    if (negative_indent_triggers.match(line.split("!", 1)[0].strip(" ")) is \
        not None):
      indent = indent[subunit_indent:]
    code[index] = "".join((indent, line))
    if (positive_indent_triggers.match(line.split("!", 1)[0].strip(" ")) is \
        not None):
      indent = "".join((subunit_indent_string, indent))

  # Separate code and comments.
  comments = ["".join(("!", line.split("!", 1)[1])) if len(line.split("!", 1)) \
      == 2 else "" for line in code]
  code = "\n".join([line.split("!", 1)[0] for line in code])

  # Put back everything between quotes.
  comment = False
  single_quotes = False
  double_quotes = False
  modified_code = ""
  for character in code:
    modified_code = "".join((modified_code, character))
    if character == "\n":
      comment = False
    elif not (single_quotes or double_quotes) and character == "!":
      comment = True
    if not (comment or double_quotes) and character == "'":
      single_quotes = not single_quotes
      if single_quotes:
        modified_code = "".join((modified_code, quotes[0]))
        del quotes[0]
    elif not (comment or single_quotes) and character == '"':
      double_quotes = not double_quotes
      if double_quotes:
        modified_code = "".join((modified_code, quotes[0]))
        del quotes[0]
  code = modified_code

  # Break lines.
  code = code.split("\n")
  for (index, line) in enumerate(code):
    length = len(line.split("\n")[- 1])
    first_linebreak = True
    while length > 80:
      left = line.split("\n")[- 1][:80]
      right = line.split("\n")[- 1][80:]
      rest = line.split("\n")[:(- 1)]
      copy = left
      if first_linebreak:
        if operator_linebreak_point.match(left):
          left = operator_linebreak_point.sub("".join((r"\1\2&\n", \
              linebreak_indent_string, r"\1&\3")), left)
        else:
          left = whitespace_linebreak_point.sub("".join((r"\1\2&\n", \
              linebreak_indent_string, r"\1&\3")), left)
        first_linebreak = False
      else:
        if operator_linebreak_point.match(left):
          left = operator_linebreak_point.sub(r"\1\2&\n\1&\3", left)
        else:
          left = whitespace_linebreak_point.sub(r"\1\2&\n\1&\3", left)
      if copy == left:
        break
      if rest:
        line = "".join(("\n".join(rest), "\n", left, right))
      else:
        line = "".join((left, right))
      length = len(line.split("\n")[- 1])
    code[index] = line

  # Put code and comments back together.
  code = "\n".join((" ".join((code[index].rstrip(" "), comments[index])) if \
      comments[index] and code[index].rstrip(" ") else "".join((code[index], \
      comments[index])) for index in range(len(code))))

  # Remove trailing whitespaces.
  code = "\n".join((entry.rstrip(" ") for entry in code.split("\n")))

  # Remove consecutive blank lines generated by unbreakable lines.
  code = re.sub(r"\n\n+", r"\n\n", code)

  # Write the output file.
  with open("/".join((output_directory, file)), "w") as output:
    output.write(code)

  # Print information.
  print(file, "is done.")
