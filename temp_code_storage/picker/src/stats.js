import React from 'react';

class Stats extends React.Component {


	render(){
		return (
		<section id="stats">
			<h2>Library stats</h2>
			<table>
				<thead>
					<tr>
						<th>Library</th>
						<th>Plate</th>
						<th>Compounds</th>
						<th>Selected <br />compounds</th>
						<th>Stat 1</th>
						<th>Stat 2</th>
						<th>Stat 3</th>
					</tr>
				</thead>
				<tbody>
					<tr>
						<td rowSpan="3">Library name 1</td>
						<td>Plate name 1</td>
						<td>666</td>
						<td>420</td>
						<td>69</td>
						<td>42</td>
						<td>21</td>
					</tr>
					<tr>
					
						<td>Plate name 2</td>
						<td>1337</td>
						<td>77</td>
						<td>13</td>
						<td>4</td>
						<td>0</td>
					</tr>
					<tr>
					
						<td>Plate name 3</td>
						<td>256</td>
						<td>128</td>
						<td>64</td>
						<td>32</td>
						<td>16</td>
					</tr>
					<tr>
						<td rowSpan="2">Library name 2</td>
						<td>Plate name 1</td>
						<td>666</td>
						<td>420</td>
						<td>69</td>
						<td>42</td>
						<td>21</td>
					</tr>
					<tr>
					
						<td>Plate name 2</td>
						<td>1337</td>
						<td>77</td>
						<td>13</td>
						<td>4</td>
						<td>0</td>
					</tr>
					<tr>
						
							<td>3 libraries</td>
							<td>5 plates</td>
							<td>4246</td>
							<td>1122</td>
							<td>228</td>
							<td>124</td>
							<td>58</td>
						
					</tr>
				</tbody>
			</table>
		</section>
		)
	}
}

export default Stats
