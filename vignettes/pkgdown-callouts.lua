-- pkgdown renders articles as minimal HTML (no Bootstrap), where Quarto
-- writes callouts as bare blockquotes and drops their type. Run this filter
-- twice: at pre-render it tags each callout with an HTML comment naming the
-- type; at post-render it wraps the blockquote carrying the tag in
-- `div.decovart-callout.callout-<type>`, which pkgdown/extra.scss styles
-- with Quarto's default appearance. Quarto still numbers each callout.
local marker_pattern = "decovart%-callout:(%a+)"

local function tag(el)
  if quarto.doc.has_bootstrap() then
    return nil
  end
  local marker = "<!-- decovart-callout:" .. (el.type or "note") .. " -->"
  local blocks = pandoc.Blocks({ pandoc.RawBlock("html", marker) })
  blocks:extend(quarto.utils.as_blocks(el.content))
  el.content = blocks
  return el
end

-- A theorem wrapped by callouty-theorem gets the generic title "Note";
-- show the numbered theorem label ("Proposition 2 (...)") there instead.
local function promote_theorem_title(blocks)
  local theorem = blocks[2]
  if not (theorem and theorem.t == "Div" and theorem.classes:includes("theorem")) then
    return
  end
  local para = theorem.content[1]
  if not (para and para.t == "Para") then
    return
  end
  local label = para.content[1]
  if not (label and label.t == "Span" and label.classes:includes("theorem-title")) then
    return
  end
  para.content:remove(1)
  while para.content[1] and (para.content[1].t == "Space" or para.content[1].t == "SoftBreak") do
    para.content:remove(1)
  end
  local rest = pandoc.utils.stringify(para.content):gsub("[%s\194\160]", "")
  if rest == "" then
    theorem.content:remove(1)
  end
  blocks[1] = pandoc.Para(label.content)
end

local function wrap(el)
  for i, block in ipairs(el.content) do
    if block.t == "RawBlock" then
      local kind = block.text:match(marker_pattern)
      if kind then
        el.content:remove(i)
        promote_theorem_title(el.content)
        return pandoc.Div({ el }, pandoc.Attr("", { "decovart-callout", "callout-" .. kind }))
      end
    end
  end
end

return {
  { Callout = tag },
  { BlockQuote = wrap }
}
