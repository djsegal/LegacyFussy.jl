function sort_lists!(main_list, other_lists...; cur_func::Function=identity)
  cur_indices = sortperm(main_list, by=cur_func)

  cur_list = main_list[cur_indices]
  _arrange_lists!(main_list, other_lists, cur_indices)
  cur_list
end

function shuffle_lists!(main_list, other_lists...)
  cur_indices = shuffle(1:length(main_list))

  cur_list = main_list[cur_indices]
  _arrange_lists!(main_list, other_lists, cur_indices)
  cur_list
end

function _arrange_lists!(main_list, other_lists, cur_indices)
  cur_lists = Any[ collect(main_list) ]
  isempty(other_lists) || append!(cur_lists, other_lists)

  for cur_list in cur_lists
      cur_list .= cur_list[cur_indices]
  end
end
