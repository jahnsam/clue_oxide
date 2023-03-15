use crate::cluster::adjacency::AdjacencyList;

/// This fuctions takes an AdjacencyList and returns a list of sets,
/// as well as a list of map indicating where vertex is.
/// Each set contains the indices of a connected subgraph.
pub fn separate_into_connected_subgraphs(adjacency_list: &AdjacencyList)
  -> ( Vec::<Vec::<usize>>, Vec::<usize> )
{

  let active_vertices = adjacency_list.get_active_vertices();
  let n_vertices = active_vertices.len();
  if n_vertices == 0{
    return (Vec::<Vec::<usize>>::new(), Vec::<usize>::new());
  }

  let mut subgraph_ids: Vec::<usize> = (0..n_vertices).collect();

  let mut counter = 1;

  for idx0 in 1..n_vertices{
    let v0 = active_vertices[idx0];

    let mut start_new_subgraph = true;
    for idx1 in (0..v0).rev(){
      let v1 = active_vertices[idx1];
      if adjacency_list.are_connected(v0,v1) {
        subgraph_ids[v0] = subgraph_ids[v1];
        start_new_subgraph = false;
        break;
      }
    }
    if start_new_subgraph{
      subgraph_ids[v0] = counter;
      counter += 1;
    }
      
  }



  let mut subgraphs: Vec::<Vec::<usize>> = (0..counter)
    .map(|_| Vec::<usize>::new()).collect();

  for (v, id) in subgraph_ids.iter().enumerate(){
    subgraphs[*id].push(v);
  }

  (subgraphs,subgraph_ids)

}

#[cfg(test)]
mod tests{
  use super::*;

  #[test]
  fn test_separate_into_connected_subgraphs(){

    // 0-1-5
    // 2-3-4
    // 6
    // 7-8
    let mut graph = AdjacencyList::with_capacity(9);
    for v in 0..9{
      graph.connect(v,v);
    }
    graph.connect(0,1);
    graph.connect(0,5);
    graph.connect(2,3);
    graph.connect(3,4);
    graph.connect(7,8);

    let expected_subgraphs = vec![
      vec![0,1,5],
      vec![2,3,4],
      vec![6],
      vec![7,8],
    ];

    let expected_ids = vec![0,0,1,1,1,0,2,3,3];

    let (subgraphs,ids) = separate_into_connected_subgraphs(&graph);

    assert_eq!(ids,expected_ids);
    assert_eq!(subgraphs,expected_subgraphs);
  }
}
