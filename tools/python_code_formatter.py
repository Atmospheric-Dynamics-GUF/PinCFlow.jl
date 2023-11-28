"""This Python code formatter unifies spacing and line breaks in Python files. \
    The spacing in comments and strings is not changed. Comments between \
    continuation lines are moved below the last of these lines. Comments \
    between quotes (i.e. in Python strings that contain line breaks) are \
    removed (a warning is given in such cases). The output does not contain \
    any extraneous whitespaces or consecutive blank lines (outside of quotes \
    and comments). The user may specify input and output directories, as well \
    as the lengths of the spacing and indent strings."""

import os
import re
import sys
import keyword

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

# Collect keywords.
keywords = re.compile(r"".join((r"\b(", r"|".join(keyword.kwlist), r")\b")))

# Collect operators.
operators = re.compile(r"(\+|\-|\*|\/|\%|\<|\>|\&|\||\~|\^|\=)")

# Define expressions for line break detection.
linebreak_comments = re.compile(r"(\n[^#\n]+\\ *(#.*)?\n)\n*(( *#.*\n)+)" \
    r"\n*([^#\n]+\\? *(#.*)?\n)")
linebreak_code = re.compile(r"(\n[^#\n]+)\\ *((#.*)?)\n([^#\n]+\\?) *" \
    r"((#.*)?\n)")

# Define expressions for line break points.
operator_linebreak_point = re.compile(r"^(\s*)(\S.*\s+)((\+|\-|\*|\/|\%|\<|\>" \
    r"|\&|\||\~|\^|\=)+\s+[^\s\+\-\*\/\%\<\>\&\|\~\^\=]+\s*)$")
whitespace_linebreak_point = re.compile(r"^(\s*)(\S.*\s+)(\S+\s*)$")
raw_single_quotes_linebreak_point = re.compile(r"^(\s*)(\S.*[^\\])(.{3})$")
raw_double_quotes_linebreak_point = re.compile(r'^(\s*)(\S.*[^\\])(.{3})$')

# Iterate over files in the input directory.
for file in os.listdir(input_directory):
  if re.search(r"\.py$", file) is None:
    continue

  # Read the input file.
  with open("/".join((input_directory, file)), "r") as input:
    code = input.read()

  # Remove trailing whitespaces.
  code = "\n".join((entry.rstrip(" ") for entry in code.split("\n")))

  # Remove consecutive blank lines.
  code = re.sub(r"\n\n+", r"\n\n", code)

  # Extract everything between quotes and raw quotes.
  comment = False
  backslash = False
  single_quotes = False
  double_quotes = False
  raw = False
  quotes = []
  raw_quotes = []
  modified_code = ""
  for (index, character) in enumerate(code):
    modified_code = "".join((modified_code, character))
    if character == "\n":
      comment = False
    elif not (single_quotes or double_quotes) and character == "#":
      comment = True
    if not (comment or backslash or double_quotes) and character == "'":
      single_quotes = not single_quotes
      if single_quotes:
        if code[index - 1] in ("r", "R"):
          raw_quotes.append("")
          raw = True
        else:
          quotes.append("")
          raw = False
    elif not (comment or backslash or single_quotes) and character == '"':
      double_quotes = not double_quotes
      if double_quotes:
        if code[index - 1] in ("r", "R"):
          raw_quotes.append("")
          raw = True
        else:
          quotes.append("")
          raw = False
    elif single_quotes or double_quotes:
      if raw:
        raw_quotes[- 1] = "".join((raw_quotes[- 1], character))
      else:
        quotes[- 1] = "".join((quotes[- 1], character))
      modified_code = modified_code[:(- 1)]
    if not backslash and character == "\\":
      backslash = True
    else:
      backslash = False
  code = modified_code

  # Join continuation lines between quotes. Comments between continuation
  # lines are removed.
  for (index, entry) in enumerate(quotes):
    if re.search(r"\n( *#.*\n)+", entry):
      entry = re.sub(r"\n( *#.*\n)+", r"\n", entry)
      print("WARNING: Comments between quotes have been removed!")
    entry = re.sub(r" *\\?\n *", r" ", entry)
    quotes[index] = entry

  # Join continuation lines between raw quotes. Comments between continuation
  # lines are removed.
  for (index, entry) in enumerate(raw_quotes):
    if re.search(r"\n( *#.*\n)+", entry):
      entry = re.sub(r"\n( *#.*\n)+", r"\n", entry)
      print("WARNING: Comments between raw quotes have been removed!")
    entry = re.sub(r"(\\?)\n", r"\1", entry)
    raw_quotes[index] = entry

  # Mark continuation lines.
  comment = False
  backslash = False
  parantheses = 0
  brackets = 0
  braces = 0
  modified_code = ""
  for character in code:
    if (not (comment or backslash) and (parantheses or brackets or braces) and \
        character == "\n"):
      character = "".join((" ", "\\", character))
    elif not comment:
      if character == "(":
        parantheses += 1
      elif character == ")":
        parantheses -= 1
      elif character == "[":
        brackets += 1
      elif character == "]":
        brackets -= 1
      elif character == "{":
        braces += 1
      elif character == "}":
        braces -= 1
      elif character == "#":
        comment = True
    elif character == "\n":
      comment = False
    if not backslash and character == "\\":
      backslash = True
    else:
      backslash = False
    modified_code = "".join((modified_code, character))
  code = modified_code

  # Join continuation lines. Comments in continuation lines are moved into the
  # first line. Comments between continuation lines are moved below the last
  # line.
  while linebreak_comments.search(code):
    code = linebreak_comments.sub(r"\1\5\3", code)
  while linebreak_code.search(code):
    code = linebreak_code.sub(r"\1\4 \2 \5", code)

  # Join implicitly concatenated raw strings. Although this seems to work,
  # non-raw strings are more complicated!
  # code = re.split(r"([rR][\"\']{2})", code)
  # count = 0
  # for index in range(3, len(code), 2):
  #   if re.fullmatch(r"\s+", code[index - 1]):
  #     code[index] = ""
  #     quote_index = (index - 1) // 2 - count
  #     raw_quotes[quote_index - 1] = "".join((raw_quotes[quote_index - 1], \
  #         raw_quotes[quote_index]))
  #     del raw_quotes[quote_index]
  #     count += 1
  # code = "".join(code)

  # Separate code and comments.
  comments = ["".join(("#", line.split("#", 1)[1])) if len(line.split("#", 1)) \
      == 2 else "" for line in code.split("\n")]
  code = "\n".join([line.split("#", 1)[0] for line in code.split("\n")])

  # Adjust whitespaces directly outside of parantheses, brackets and braces.
  code = re.sub(r"\) *", ") ", code)
  code = re.sub(r"(\S) *\(", r"\1(", code)
  code = re.sub(r"\] *", "] ", code)
  code = re.sub(r"(\S) *\[", r"\1[", code)
  code = re.sub(r"\} *", "} ", code)
  code = re.sub(r"(\S) *\{", r"\1{", code)

  # Add whitespaces behind keywords.
  code = keywords.sub(r"\1 ", code)

  # Add whitespaces around single operators.
  code = operators.sub(r" \1 ", code)

  # Adjust whitespaces in power notation.
  code = re.sub(r"\b([\d.]+) *(e) *((\+|\-)?) *([\d.]+)\b", r"\1\2\3\5", code, \
      flags = re.IGNORECASE)

  # Adjust whitespaces in decorators.
  code = re.sub(r"(\n *@) *", r"\1", code)

  # Adjust whitespaces at combined operators.
  code = re.sub(r"\/ +\/", "//", code)
  code = re.sub(r"\* +\*", "**", code)
  code = re.sub(r"\= +\=", "==", code)
  code = re.sub(r"\! +\=", "!=", code)
  code = re.sub(r"\< +\=", "<=", code)
  code = re.sub(r"\> +\=", ">=", code)
  code = re.sub(r"\< +\<", "<<", code)
  code = re.sub(r"\> +\>", ">>", code)
  code = re.sub(r"\+ +\=", "+=", code)
  code = re.sub(r"\- +\=", "-=", code)
  code = re.sub(r"\* +\=", "*=", code)
  code = re.sub(r"\/ +\=", "/=", code)
  code = re.sub(r"\% +\=", "%=", code)
  code = re.sub(r"\& +\=", "&=", code)
  code = re.sub(r"\| +\=", "|=", code)
  code = re.sub(r"\^ +\=", "^=", code)
  code = re.sub(r"\- +\>", "->", code)
  code = re.sub(r"\< +\-", "<-", code)

  # Adjust whitespaces around colons.
  code = ":".join(("".join((" ", entry.strip(" "))) if index and (not \
      ":".join(code.split(":")[:index]).split("{")[- 1].count("}") and \
      ":".join(code.split(":")[:index]).split("{")[- 1].count("[") \
      == ":".join(code.split(":")[:index]).split("{")[- 1].count("]")) else \
      entry.strip(" ") for (index, entry) in enumerate(code.split(":"))))

  # Adjust whitespaces around dots, commas and semicolons.
  code = re.sub(r" *\. *", ".", code)
  code = re.sub(r" *, *", ", ", code)
  code = re.sub(r" *; *", "; ", code)

  # Remove whitespaces directly inside parantheses.
  code = re.sub(r" +\)", ")", code)
  code = re.sub(r"\( +", "(", code)
  code = re.sub(r" +\]", "]", code)
  code = re.sub(r"\[ +", "[", code)
  code = re.sub(r" +\}", "}", code)
  code = re.sub(r"\{ +", "{", code)

  # Remove extraneous whitespaces.
  code = re.sub(r"(\S) +", r"\1 ", code)

  # Adjust subunit indents.
  if re.search(r"\n(class|def|elif|else|except|finally|for|if|try|while|with)" \
      r".*\n+( *)\S", code):
    match = re.search(r"\n(class|def|elif|else|except|finally|for|if|try|" \
        r"while|with).*\n+( *)\S", code)
    code = re.sub(match.expand(r"\2"), subunit_indent_string, code)

  # Put back everything between quotes and raw quotes.
  comment = False
  single_quotes = False
  double_quotes = False
  modified_code = ""
  for (index, character) in enumerate(code):
    modified_code = "".join((modified_code, character))
    if character == "\n":
      comment = False
    elif not (single_quotes or double_quotes) and character == "#":
      comment = True
    if not (comment or double_quotes) and character == "'":
      single_quotes = not single_quotes
      if single_quotes:
        if code[index - 1] in ("r", "R"):
          modified_code = "".join((modified_code, raw_quotes[0]))
          del raw_quotes[0]
        else:
          modified_code = "".join((modified_code, quotes[0]))
          del quotes[0]
    elif not (comment or single_quotes) and character == '"':
      double_quotes = not double_quotes
      if double_quotes:
        if code[index - 1] in ("r", "R"):
          modified_code = "".join((modified_code, raw_quotes[0]))
          del raw_quotes[0]
        else:
          modified_code = "".join((modified_code, quotes[0]))
          del quotes[0]
  code = modified_code

  # Break lines.
  code = code.split("\n")
  for (index, line) in enumerate(code):
    length = len(line.split("\n")[- 1])
    first_linebreak = True
    backslash = False
    single_quotes = False
    double_quotes = False
    raw = False
    while length > 80:
      left = line.split("\n")[- 1][:80]
      right = line.split("\n")[- 1][80:]
      rest = line.split("\n")[:(- 1)]
      copy = left
      # Check if the left part ends between quotes.
      for (position, character) in enumerate(left[:(- 3)]):
        if not (backslash or double_quotes) and character == "'":
          single_quotes = not single_quotes
          if single_quotes:
            if left[position - 1] in ("r", "R"):
              raw = True
            else:
              raw = False
          else:
            raw = False
        elif not (backslash or single_quotes) and character == '"':
          double_quotes = not double_quotes
          if double_quotes:
            if left[position - 1] in ("r", "R"):
              raw = True
            else:
              raw = False
          else:
            raw = False
        if not backslash and character == "\\":
          backslash = True
        else:
          backslash = False
      # Insert the first linebreak.
      if first_linebreak:
        if single_quotes and raw:
          left = raw_single_quotes_linebreak_point.sub("".join((r"\1\2' \\\n", \
              linebreak_indent_string, r"\1r'\3")), left)
          single_quotes = False
        elif double_quotes and raw:
          left = raw_double_quotes_linebreak_point.sub("".join((r'\1\2" \\\n', \
              linebreak_indent_string, r'\1r"\3')), left)
          double_quotes = False
        else:
          if operator_linebreak_point.match(left):
            left = operator_linebreak_point.sub("".join((r"\1\2\\\n", \
                linebreak_indent_string, r"\1\3")), left)
          else:
            left = whitespace_linebreak_point.sub("".join((r"\1\2\\\n", \
                linebreak_indent_string, r"\1\3")), left)
        first_linebreak = False
      # Insert subsequent linebreaks.
      else:
        if single_quotes and raw:
          left = raw_single_quotes_linebreak_point.sub(r"\1\2' \\\n\1r'\3", \
              left)
          single_quotes = False
        elif double_quotes and raw:
          left = raw_double_quotes_linebreak_point.sub(r'\1\2" \\\n\1r"\3', \
              left)
          double_quotes = False
        else:
          if operator_linebreak_point.match(left):
            left = operator_linebreak_point.sub(r"\1\2\\\n\1\3", left)
          else:
            left = whitespace_linebreak_point.sub(r"\1\2\\\n\1\3", left)
      # Update the split.
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