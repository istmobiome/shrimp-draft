local collected_links = {}
local seen_links = {}

-- Collect unique links
local function add_link(link_elem)
  local link_key = pandoc.utils.stringify(link_elem.content) .. "|" .. link_elem.target
  if not seen_links[link_key] then
    seen_links[link_key] = true
    table.insert(collected_links, link_elem:clone())
  end
end

-- First pass: collect links
local function collect_links(el)
  if el.t == "Link" then
    add_link(el)
  end
  return nil
end

-- Second pass: append a Collected Links section at the end
local function add_links_section(doc)
  if #collected_links == 0 then
    return doc
  end

  local links_header = pandoc.Header(4, "Collected Links", pandoc.Attr("", {"appendix"}))
  local links_list_items = {}
  for _, link in ipairs(collected_links) do
    table.insert(links_list_items, pandoc.Plain({ link }))
  end
  local bullet_list = pandoc.BulletList(links_list_items)

  table.insert(doc.blocks, links_header)
  table.insert(doc.blocks, bullet_list)

  return doc
end

return {
  { Link = collect_links },
  { Pandoc = add_links_section }
}
