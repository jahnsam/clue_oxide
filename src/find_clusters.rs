use crate::signal;

pub struct Cluster{
  vertices: Vec::<usize>,
  //subclusters Option<&Vec::<Cluster> >, // Maybe?
  signal: Option<Signal>,
  auxiliary_signal: Option<Signal>,
}

impl ToString for cluster{
  fn to_string(&self) -> String {
    // "[1,2,3,4]"
    to_cluster_string(self.vertices)
  }
}
pub fn to_cluster_string(cluster: &[usize]){

}
/*
pub struct Clusters{
  clusters: Vec::<HashMap::<Cluster>>,
  //subcluster_info // Maybe?
}
*/

pub fn find_clusters(
    andjacency_matrix: Matrix, // what data structure?
    max_size: usize) 
  -> Result< Vec::<HashMap::<Cluster>>, CluEError>
{

  let mut clusters = Vec::<HashMap<Cluster>>::with_capacity(max_size);

  let one_clusters = HashMap::new(); 
  // Identify all 1-clusters, as those with adjacency_matrix[[ii,ii]] == true.
  for ii in 0 .. adjacency_matrix.n_rows(){
    if !adjacency_matrix[[ii,ii]] {continue}

    let cluster = Cluster::new(vec![ii]));
    let key cluster.to_string();
    one_clusters.entry(key).or_insert(cluster);
  }


  clusters.push(one_clusters)

  for isize in 1..=max_size{

    // Build n-clusters from (n-1)-clusters.
    let n_clusters = build_n_clusters(&mut clusters[isize -1],adjacency_matrix);


  }
}

//------------------------------------------------------------------------------

fn build_n_clusters(
    &mut n_minus_1_clusters: &[Cluster], adjacency_matrix: &Matrix)
  -> Result<HashMap::<Cluster>,CluEError>
{

  let mut new_clusters = HashMap::new();

  //https://stackoverflow.com/questions/45724517/how-to-iterate-through-a-hashmap-print-the-key-value-and-remove-the-value-in-ru
  for (key, value) in &*n_minus_1_clusters {

    let neighbors = find_neighbors(value,adjacency_matrix);

    for vertex in neighbors{
    
      let new_cluster = sort(concatenate(value,vertex));

      let key new_cluster.to_string();
      new_clusters.entry(key).or_insert(new_cluster)
    } 

  }

  Ok(new_clusters)
}

