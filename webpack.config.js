module.exports = {
    module: {
        rules: [
            {
                test: /\.js$/,
                exclude: /node_modules/,
                use: {
                    loader: "babel-loader"
                }
            },
            {
				 test: /\.(jpg|png|svg|gif)$/,
				 loader: 'file-loader',
				 options: {
					 name: '[path][name].[ext]'
				 }
			},
 
        ]
    }
}
