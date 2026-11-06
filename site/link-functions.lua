-- Links inline code naming an exported function, e.g. `spatial_similarity()`,
-- to its section of the Reference
local functions = { spatial_similarity = true, plot_neighborhoods = true }

return {
  {
    traverse = "topdown",
    -- Leave headings and existing links alone
    Header = function(el) return el, false end,
    Link = function(el) return el, false end,
    Code = function(el)
      local name = el.text:match("^([%w_]+)%(%)$")
      if name and functions[name] then
        return pandoc.Link(el, "reference.qmd#" .. name), false
      end
    end,
  },
}
