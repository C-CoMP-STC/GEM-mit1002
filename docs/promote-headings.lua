-- PDF build only (docs/build.sh).
--
-- continuous-curation.md keeps its own top-level "# Continuous Curation..."
-- heading so the file reads properly on GitHub. In the PDF the title comes
-- from metadata.yaml instead, so that heading would be a duplicate -- and
-- worse, with --number-sections it becomes section 1 and every real section
-- (Motivation, Set-Up, The Continuous Curation Loop...) becomes a subsection
-- of it.
--
-- So: drop the first level-1 heading, and promote everything else one level,
-- making the "##" sections top-level in the PDF.
--
-- Delete --lua-filter=promote-headings.lua from build.sh if you ever want the
-- document's own title heading to appear in the PDF as well.

local dropped = false

function Header(el)
  if not dropped and el.level == 1 then
    dropped = true
    return {}
  end
  el.level = math.max(1, el.level - 1)
  return el
end
